# ortho2tree

## SETTING UP AN ANALYSIS

This is done via the creation of a YAML configuration file, copying the example `test.cfg` provided, renaming it to a meaningful name, which will be the `set` name.

For example, to set up an analysis on a number of `plants` proteomes, we could:
```cp test.cfg plants.cfg```

and then edit the config file `plants.cfg` to specify the organisms to analyze.

This is done by editing the definition of the `tax2org` dictionary which could become, in our example:
```
tax2org:
    3702: arabidopsis
    4577: maize
    39947: rice
    3635: cotton
    3880: barrel_medic
    4097: tobacco
    4565: wheat
    112509: barley
```

**Note** -- the names listed after the taxon identifier need to be the names used by [Panther](https://www.pantherdb.org/) database. They are used to to fetch the Panther orthogroup accession mapping.

## RUNNING THE ANALYSIS

The analysis is started by launching the command (in our example with the plants `set`):
```./ortho2tree.py -set plants```

The pipeline will first load [Panther](https://www.pantherdb.org/) mapping data from the [Panther ftpsite](ftp://ftp.pantherdb.org/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/) and [UniProtKB](https://www.uniprot.org) Gene-Centric data for the organisms in the analysis. Both sources of data are cached to local files.

In addition to the Panther ortholog mapping, the script needs the
mapping of isoforms to canonicals (Panther only provides canonical
accessions); this is provided by Gene-Centric. Gene-Centric data is
either retrieved via uniprot database access (internal use), or it can
be provided as a local or web based tsv file (available on request,
please contact us) or alternatively it can be constructed by parsing
the `.gene2acc` files for the Reference Proteomes, [available from
ftp](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/).

These two sources of data are processed and combined into a single
python [pandas](https://pandas.pydata.org/) DataFrame, which can be cached to a file (config file option `use_cached_orthogroup_data`). That file can later be used for subsequent runs (e.g. trying different parameters), controlled by two other config file options (`dump_orthogroup_data` and `instamp`).

The pipeline will work (optionally in parallel, depending on the
number of `threads` specified in the config file) on all the
orthogroups identified, after discarding those with a number of
proteomes less than the specified `min_taxa_threshold` (option in config).

For each orthogroup, it will:
- fetch sequences (optionally saving them under the directory specified as `fasta_dir` in config)
- create a multiple sequence alignment (optionally saving it under the directory specified as `aln_data_dir`)
- build a gap-based Neighbour-Joining tree (optionally saving it under `tree_data_dir`)
- identify and optionally save to disk (under `n_data_dir`) low-cost clades
- rank the best low-cost clades to confirm existing canonicals or suggest replacements (the output of the pipeline)

## OUTPUT FILES

If the pipeline is run on a single group of orthologs (or a list of
few groups) with a command line option like `-id PTHR12345:SF67`, the
output will be printed.

If the `-id orthogroup` option is not specified, all the Panther
orthogroups available for those organisms will be analysed, and the
results for each orthogroup will be written to four files (timestamped
according to the `outstamp` parameter in the config file or with the -outstamp command line option):

- `output_changes`: the proposed canonical replacements
- `output_confirm`: the clades with confirmed canonicals
- `output_skipped`: groups which did not propose changes nor included confirmed clades, usually due to now good low-cost clades being found
- `output_gc`: similar to `output_changes` but with md5 checksums for sequences and excluding previously made suggestions listed a `prevgc_file` (if specified in the config)

Each of the `output_` files has a header (lines starting with #) with the parameters used for the analysis. The last line of the header has the column names and begins with the text `# pthr_id`.

## FORMAT OF `output_changesyymmdd` (`output_changes` + `yymmddd` timestamp)

Each proposed change is in a line with the following tab-separated fields:

```
pthr_id           Panther orthogroup identifier
taxon             organism identifier as (OSCODE)[https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/docs/speclist.txt]
canon_acc         UniProtKB accession of the current canonical (sp: reviewed, tr: unreviewed)
canon_len         sequence length of the current canonical
prop_acc          UniProtKB accession of the proposed change (sp: reviewed, tr: unreviewed)
prop_len          sequence length of the proposed change
rank_score        score weighted by prop_cost, n_sp, w_canon, etc. used to rank proposed changes
canon_cost        cost of clade with the specified taxa using the canonical sequences from those taxa
prop_cost         cost of clade if proposed canonicals replace canonicals
n_sp              number of reviewed (UniProtKB/SwissProt) sequences
n_tax             number of taxa in the clade
n_canon           number of canonical sequences in clade
wn_canon          weighted number of canonical sequences in clade (based on taxon preference: human and mouse count 2X, rat 1X, other canonicals are not counted)
scaled_prop_f     0-1 exp scaled proportion of isoforms in clade
scaled_p_cost     0-1 exp scaled prop_cost
p_len_diff        0-1 exp scaled difference in length
clade             name of the clade from the Neighbour-Joining tree
clade_members     all the members of the clade (comma separated list of OSCODE:ACCESSION identifiers)
MANEstatus        agreement of clade with MANE (if human in analysis): MANE_good, MANE_bad, NAM (orthogroup clade has no HUMAN sequence, so no MANE assignment)
```

## FORMAT OF `output_confirmyymmdd` (`output_confirm` + `yymmdd` timestamp):

```
pthr_id           Panther orthogroup identifier
canon_cost        cost of the clade
n_sp              number of reviewed (UniProtKB/SwissProt) sequences
n_tax             number of taxa in the clade
clade             name of the clade from the Neighbour-Joining tree
clade_members     all the members of the clade (comma separated list of OSCODE:ACCESSION identifiers)
MANEstatus        agreement of clade with MANE (if human in analysis): MANE_good, MANE_bad, NAM (orthogroup clade has no HUMAN sequence, so no MANE assignment)
```

## FORMAT OF `output_skippedyymmdd` (`output_changes` + `yymmdd` timestamp)

```
pthr_id           Panther orthogroup identifier
clade             name of the clade from the Neighbour-Joining tree
reason            reason for the clade being skipped
info              additional information (for example the cost of the clade)
```

Most orthogroups will either be in `output_changes` or
`output_confirm`. Some orthogroups may be in both `output_changes`
and `output_confirm` if an orthogroup produces multiple clades with
both proposed changes and confirmed canonicals. Clades are listed in
`output_skipped` only if there are no other clades in the orthogroup
with either confirmed or proposed changes.
