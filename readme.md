# DMS Primer Design

Code to generate genetic variants of a gene of interest.

The `primer_design.py` requires a `config.yaml` and `.gb` file for a gene of interest to generate a list of doped primers.  The designed primers can be used to produce a library of single-residue variants that can be stained by a single nucleotide change (or SNP) in the gene of interest.

## Quick Start
```
conda env create -f environment.yaml
conda activate dms_primer_design_env
python 
```
