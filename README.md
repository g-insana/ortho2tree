# ortho2tree

The [UniProt](https://www.uniprot.org) [Reference Proteomes dataset](https://www.uniprot.org/help/reference_proteome) seeks to provide complete proteomes for an evolutionarily diverse, less redundant, set of organisms. 

As higher eukaryotes often encode multiple isoforms of a protein from a single gene, the Reference Proteome pipeline selects a single representative (‘canonical’) sequence. UniProt identifies canonical isoforms using a ‘Gene-Centric’ approach: proteins are grouped by gene-identifier and for each gene a single protein sequence is chosen. 

For unreviewed (UniProtKB/TrEMBL) protein sequences (and for some reviewed sequences), the longest sequence in the Gene-Centric group is usually chosen as canonical. This can create inconsistencies, selecting canonical sequences with dramatically different lengths for orthologous genes.

The Ortho2tree data pipeline examines Gene-Centric canonical and isoform sequences from sets of orthologous proteins (from PantherDB), builds multiple alignments, constructs gap-distance trees, and identifies low-cost clades of isoforms with similar lengths. Canonical choices can be either confirmed or a better one proposed.

An overview of the pipeline is shown in this figure:
![ortho2tree pipeline overview](ortho2tree_pipeline.jpg)

The pipeline can retrieve protein sequences using direct access to the UniProt databases or using the UniProt web API.

Data processing is done on DataFrames employing vectorization and all the orthogroups are processed in parallel:
- building Multiple Sequence Alignments (via muscle)
- calculating gap-based Neighbour-Joining trees (via BioPython with a modified pairwise distance function)
- scanning trees to identify low-cost clades
- ranking best low-cost clades to confirm existing canonicals or suggest replacements

# Contents of the repository
```
ortho2tree.py    # main script to use to run the pipeline on the command line
ortho2tree.ipynb # jupyter notebook to run the pipeline interactively
ortho2tree/      # modules folder
requirements.txt # list of needed packages
README.md        # this text
test/            # folder containing data ready for a test run
test.cfg         # configuration file for a test run
qfomam           # data for the analsys of 2022_05 data on QfO mammals
```

## INSTALLATION
- git clone the repository: 

```git clone https://github.com/g-insana/ortho2tree.git``` 
- install requirements (via conda/mamba or via pip): 

```cd ortho2tree && pip install -r requirements.txt # example installation via pip``` 

## TEST USAGE
- test run of a single group
```./ortho2tree.py -set test -id PTRH43715:SF1```
- full analysis run of a set
```./ortho2tree.py -set test -id PTRH43715:SF1```
