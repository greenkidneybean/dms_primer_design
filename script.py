#!/usr/bin/env python3

# imports
import argparse
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio import SeqIO
import math
import numpy as np
import pandas as pd
import main_package # my package

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("wt", help="Genbank file path containing wild type (WT) sequence", type=str)
parser.add_argument("o", help="Output prefix", type=str)
parser.add_argument("--vector", help="Genbank file path containing vector sequence", type=str, default=False)
parser.add_argument("--codon_table", help="Specify codon table to use", type=str, default='Standard')
parser.add_argument("--homo_len", help="Length of homology arm in fwd primer", type=int, default=20)
parser.add_argument("--oligo_len", help="Ideal max total length of oligo", type=int, default=60)
parser.add_argument("--melt_temp", help="Melting temp of fwd primer", type=int, default=50)
parser.add_argument("--rev_melt_temp", help="Melting temp of rev primer", type=int, default=55)
parser.add_argument("--syn_snp_rate", help="Percentage of synonymous SNPs 0-1", type=float, default=.05)
parser.add_argument("--remove_stop_rate", help="Percentage of stop codon SNPs to remove, default = 90% of stop SNPs", type=float, default=.90)
parser.add_argument("--rng_seed", help="Set seed for repoducibly selecting synonymous codon sites", type=int, default=42)
parser.add_argument("--out_dir", help='Local output directory e.g. "data"', type=str)
args = parser.parse_args()

# parse genbank files
wt_file = SeqIO.read(args.wt, 'genbank')

# check for vector file
if not args.vector:
    args.vector = args.wt
vector_file = SeqIO.read(args.vector, 'genbank')

wt_seq = str(wt_file.seq.upper())
vector_seq = str(vector_file.seq.upper())

# ERROR CHECKS
if len(wt_seq) != len(vector_seq):
    print('ERROR: WildType and Vector GenBank sequences are not of equal length')
    exit()
# check for -20 bp homology
# check that the strand is going forward

# get start and stop of gene for codon positions
for feature in wt_file.features:
    if feature.type == 'gene':
        gene_start = feature.location.start.position
        gene_end = feature.location.end.position

# setup seq_data
seq_data = {}
seq_data['wt_seq'] = wt_seq
seq_data['vector_seq'] = vector_seq
seq_data['gene_start'] = gene_start
seq_data['gene_end'] = gene_end
seq_data['fasta_file'] = []
seq_data['df'] = pd.DataFrame()
seq_data['rng'] = np.random.RandomState(int(args.rng_seed))

# this needs to be fixed (user input? yaml?)
targ_windows = ['window_1', 'window_2', 'window_3']

for feature in wt_file.features:
    if feature.type not in targ_windows:
        continue

    start_index = feature.location.start.position
    window_end = feature.location.end.position

    # loop for each sub_window
    sub_window_n = 1
    while start_index < window_end: # this could be an issue to toggle
        data_dict = {}
        data_dict['start_index'] = start_index
        data_dict['window_end'] = window_end
        data_dict['sub_window_name'] = f"{str(feature.type)}-{sub_window_n}"

        # 1. homology arm
        data_dict = main_package.primer_design.homology_arm(seq_data, data_dict, args)

        # 2. reverse primer
        data_dict = main_package.primer_design.reverse_primer(seq_data, data_dict, args)

        # 3. forward primer
        data_dict = main_package.primer_design.forward_primer(seq_data, data_dict, args)

        # 4. variant window
        seq_data, data_dict = main_package.primer_design.sub_window(seq_data, data_dict, args)

        # reset the start index for the next mini-window
        start_index = data_dict['primer_start']
        sub_window_n += 1

# setup .fa output, truncate if file exists
file = open(f"{args.o}.fa",'w+')
file.writelines(seq_data['fasta_file'])
file.close()

# polish dataframe
df = seq_data['df']
df['position'] = df['position'].astype(int)

df['forward_primer_tm'] = df['forward_primer'].apply(lambda x: mt.Tm_NN(x)).round(1)
df['forward_primer_gc'] = df['forward_primer'].apply(GC).round(1)
df['forward_primer_len'] = df['forward_primer'].str.len()

df['reverse_primer_tm'] = df['reverse_primer'].apply(lambda x: mt.Tm_NN(x)).round(1)
df['reverse_primer_gc'] = df['reverse_primer'].apply(GC).round(1)
df['reverse_primer_len'] = df['reverse_primer'].str.len()

cols = ['name','sub_window_name','wt_codon','position','iupac_codon','codon_sub','iupac_aa','add_synonymous_codon','contains_missense_stop','remove_missense_stop_codon','primer','homology_arm','sub_window','forward_primer','forward_primer_tm','forward_primer_gc','forward_primer_len','reverse_primer_name','reverse_primer','reverse_primer_tm','reverse_primer_gc','reverse_primer_len']
df = df[cols]

# save dataframe as .tsv
df.to_csv(f'{args.o}.tsv', index=False, sep='\t')
