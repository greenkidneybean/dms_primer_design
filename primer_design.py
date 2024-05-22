#!/usr/bin/env python3

# need to correct the output_path if it doesn't exist

# imports
import os
import argparse
import yaml
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC
import math
import numpy as np
import pandas as pd
import dms_primer_design # custom package

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("config", help="configuration file", type=str)
parser.add_argument("--default_config", help="custom default configuration file", type=str, default=None)
args = parser.parse_args()

# parse and merge config files
cfg = dms_primer_design.primer_design.merge_config_files(args.config, args.default_config)

print()
print("CONFIGURATION:")
for key, value in cfg.items():
    print(f"{key}: {value}")
print()

# parse genbank files
wt_file_path = cfg['wt_seq']
wt_file = SeqIO.read(wt_file_path, 'genbank')

# check for vector file
vector_file_path = cfg['vector_seq']
if not vector_file_path:
    vector_file_path = wt_file_path
vector_file = SeqIO.read(vector_file_path, 'genbank')

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
# this should be a class
seq_data = {}
seq_data['wt_seq'] = wt_seq
seq_data['vector_seq'] = vector_seq
seq_data['gene_start'] = gene_start
seq_data['gene_end'] = gene_end
seq_data['fasta_file'] = []
seq_data['df'] = pd.DataFrame()
seq_data['rng'] = np.random.RandomState(cfg['rng_seed'])

# identify gene "windows" where variants will be made
targ_windows = cfg['variant_windows']

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
        data_dict = dms_primer_design.primer_design.homology_arm(seq_data, data_dict, cfg)

        # 2. reverse primer
        data_dict = dms_primer_design.primer_design.reverse_primer(seq_data, data_dict, cfg)

        # 3. forward primer
        data_dict = dms_primer_design.primer_design.forward_primer(seq_data, data_dict, cfg)

        # 4. variant window
        seq_data, data_dict = dms_primer_design.primer_design.sub_window(seq_data, data_dict, cfg)

        # reset the start index for the next mini-window
        start_index = data_dict['primer_start']
        sub_window_n += 1

# setup .fa output, truncate if file exists
output_prefix = cfg['output_prefix']
output_path = cfg['output_path']

# create output path
if not os.path.exists(output_path):
    os.makedirs(output_path)

file = open(f"{output_path}{output_prefix}.fa",'w+')
file.writelines(seq_data['fasta_file'])
file.close()

# polish dataframe
df = seq_data['df']
df['position'] = df['position'].astype(int)

df['primer_len'] = df['primer'].str.len()
df['forward_primer_tm'] = df['forward_primer'].apply(lambda x: mt.Tm_NN(x)).round(1)
df['forward_primer_gc'] = df['forward_primer'].apply(GC).round(1)
df['forward_primer_len'] = df['forward_primer'].str.len()

df['reverse_primer_tm'] = df['reverse_primer'].apply(lambda x: mt.Tm_NN(x)).round(1)
df['reverse_primer_gc'] = df['reverse_primer'].apply(GC).round(1)
df['reverse_primer_len'] = df['reverse_primer'].str.len()
df['wt_aa'] = df['wt_codon'].apply(lambda x: str(Seq(x).translate()))

cols = ['name','sub_window_name','wt_codon','wt_aa','position','iupac_codon','iupac_aa','iupac_sub','codon_subs','add_synonymous_codon','contains_missense_stop','remove_missense_stop_codon','primer_len','primer','homology_arm','sub_window','forward_primer','forward_primer_tm','forward_primer_gc','forward_primer_len','reverse_primer_name','reverse_primer','reverse_primer_tm','reverse_primer_gc','reverse_primer_len']
df = df[cols]

total_variants = df.iupac_aa.str.len().sum()
syn_variants = df.query('add_synonymous_codon == 1').position.nunique()
syn_location = df.query('add_synonymous_codon == 1').position.unique().tolist()
miss_variants = total_variants - syn_variants
total_stops = df.query('contains_missense_stop == 1').position.nunique()
stops_made = df.query('contains_missense_stop == 1').query('remove_missense_stop_codon == 0').position.nunique()
stop_location = df.query('contains_missense_stop == 1').query('remove_missense_stop_codon == 0').position.unique().tolist()

# print report
print('Total variants:', total_variants)
print('Missense variants:', miss_variants)
print('Synonymous variants:', syn_variants)
print('Synonymous variant location:', syn_location)
print('Fraction of stop variants:', f'{stops_made}/{total_stops}')
print('Stop variant location:', stop_location)
print('Doped variant primers:', df.shape[0])
print('Variant sub-windows:', df.sub_window_name.nunique())

# save dataframe as .tsv
df.to_csv(f'{output_path}{output_prefix}.tsv', index=False, sep='\t')
