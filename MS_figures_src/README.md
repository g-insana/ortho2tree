## **Ortho2tree** manuscript figure/plot files

(3-March-2024)

This archive contains all of the datafiles and .R code to recreate the figures in the manuscript: 
*Improved selection of canonical proteins for reference proteomes*.

To recreate pdf files for the figures, run: `plot_MS_pub.sh`

Each `.R` script has an associated figure number associated with
it. `fig1_percid_gaps3.R` is used for Fig. 1 and Suppl. Fig. 4.  In
addition, several of the scripts use `set_color_a90.R` to set line
types, symbols, and colors.

In addition to the .R scripts to create the files, various data files are included:

- `aln_data/` contains `*.summ2i` and `*.summ2i_m` files that report the
best hits for various search runs.  For example,
`human_vs_mam_canon.summ2i` contains the best hits for the searches of
human canonical proteins vs gorgo (gorilla), mouse, rat, bovin, and
mondo (opossum) reference proteomes. `*.summ2i` files are plotted by
`fig1_percid_gaps3.R` for Fig. 1 and suppl. Fig. 4.
`aln_data/*.summ2i_m` are derived from `*.summ2i` files, but have an
an additional file that indicates the Panther orthogroup for the human
(or mouse) query.  These files are used by `fig5_percid_gaps6.R`.

- `cnt_data/` contains a set of `.tab` files that report the clade sizes
for each of the orthogroups for a specific UniProtKB release, used in
Fig. 7.  Thus, `mam2023_01_clade_stats2.tab` contains the clade sizes
for the proposed changes made starting with the UP2023_01 release
set. The Trembl/Trembl proposed changes from this analysis were
actually incorporated into the UP2023_02 release.  The ortho2tree
changes were proposed based on the 8 QfO mammals in release UP2022_05;
the `qfomam_2022_05_clade_stats2.tab` file has the clade size
statistics for that analysis (which was actually incorporated into UP
release 2023_01, the first release using ortho2tree suggestions).

*Note that these folders are compressed as .tar.gz archives, and will be uncompressed when running the script `plot_MS_pub.sh`*

Other data files include:

- `qfomam_output_changes240206.one` and `qfomam_output_confirm240206`,
  used to create Suppl. Fig. 1

- `pmam2023_05_output_gc211004` used to create Suppl. Fig 3

- `history_new_changes_by_species.tsv` used to create Fig. 6

- `*.times` -- files that provide the divergence times used in Figs. 1 and Suppl. Fig. 4
