# DMS Primer Design

Code that designes primers used to generate genetic variants from a gene of interest.

The `primer_design.py` script requires a `config.yaml` and `.gb` file for a gene of interest, and outputs a `.tsv` file of doped primers.  This script specifically designs primers that generate single-residue missense vairants that are attained by a single nucleotide change (or SNP) from the original gene of interest sequence. Example inputs are provided in the `test_data` directory. Variants are generated across "tiles" of the coding sequence that are listed as Features in the `.gb` and listed under "variant_windows" in the `config.yaml`.  Note that the variant codons are flanked with synonymous changes to the "vector_seq" file.

## Quick Start
```
conda env create -f environment.yaml
conda activate dms_primer_design_env
python 
```

## Config File
output_prefix: the name of the gene of interest 
output_path: directory where `.tsv` file will be saved 
rng_seed: integer used to seed where synonymous and stop variants are designed 
wt_seq: path to `.gb` file containing the coding sequence for the gene of interest from which SNP-accesible variants are designed 
vector_seq: path to alternate `.gb` file containing the coding sequence for the gene of interest.  Can be identical to the wt_seq path, but provides the option to design variants for a wildtype nucleotide sequence that may differ from a codon-optimized nucleotide sequence on the vector plasmid. 
homology_length: length of homology arms for Gibson Assembly of PCR fragments 
max_oligo_length: the ideal maximum length for each primer, but will extend if necessary to increase the melting temperature of the primer 
fwd_primer_melt_temp: melting temperature of the forward doped primer that generates the variant 
rev_primer_melt_temp: melting temperature of the reverse primer for a given tile 
mutagenesis_type: default is "missense" to generate only SNP-accessible variants 
variant_windows: a list of feature names in the `.gb` file that specify where each tile of variants should be made 
synonymous_variant_rate: float value ranging from 0-1 denoting the percent of synonymous variants, e.g a value of .05 would result in 5% of the designed primers encoding synonymous variants, used to generate control variants 
remove_stop_variant_rate: float value ranging from 0-1 denoting the percent of stop variants removed from primer design, e.g a value of 1.0 would result in all stop variants being excluded from the primer design, used to generate control variants 
