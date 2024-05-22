#!/usr/bin/env python3

import argparse
import pandas as pd

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("primer_files", help=".tsv primer files to be sorted", type=str, nargs='+')
parser.add_argument("--out", help="output file name", type=str, default="sorted_primers.tsv")
parser.add_argument("--seperator", help="integer value to separate primer lengt", type=int, default=60)
args = parser.parse_args()

# merge primer input files
df_list = [pd.read_csv(primer_df, sep='\t') for primer_df in args.primer_files]

df = pd.concat(df_list)
# primers less then or equal to sep_length value
df1 = df.query(f'primer_len <= {args.seperator}')
df2 = df.query(f'primer_len > {args.seperator}')

df = pd.concat([df1,df2])
df.to_csv(args.out, index=False, sep='\t')
