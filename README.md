# ortho2tree

The [UniProt](https://www.uniprot.org) [Reference Proteomes dataset](https://www.uniprot.org/help/reference_proteome) seeks to provide complete proteomes for an evolutionarily diverse, less redundant, set of organisms. 

As higher eukaryotes often encode multiple isoforms of a protein from a single gene, the Reference Proteome pipeline selects a single representative (‘canonical’) sequence. UniProt identifies canonical isoforms using a ‘canonicalGene-Centric’ approach: proteins are grouped by gene-identifier and for each gene a single protein sequence is chosen. 

For unreviewed (UniProtKB/TrEMBL) protein sequences (and for some reviewed sequences), the longest sequence in the Gene-Centric group is usually chosen as canonical. This can create inconsistencies, selecting sequences with dramatically different lengths as canonical for orthologous genes. Biologically, it is unlikely that orthologous mammalian proteins differ greatly in length. 

The Ortho2tree data pipeline examines Gene-Centric canonical and isoform sequences from sets of orthologous proteins (from PantherDB), and suggests replacements for canonicals that have lengths very different from closely related orthologs. 

# Contents of the repository
```
ortho2tree.py    # main script to use to run the pipeline on the command line
ortho2tree.ipynb # jupyter notebook to run the pipeline interactively
requirements.txt # list of needed packages
config_muscle.py # to configure path to muscle executable and muscle version used.
o2t_*.py         # source code modules used by the ortho2tree pipeline
README.md        # this text
test/            # folder containing data ready for a test run
test.cfg         # configuration file for a test run
```

## INSTALLATION
- git clone the repository: 

```git clone https://github.com/g-insana/ortho2tree.git``` 
- install requirements (via conda/mamba or via pip): 

```cd ortho2tree && pip install -r requirements.txt # example installation via pip``` 

## USAGE
- test run of a single group 
- full analysis run of a set 
