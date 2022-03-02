# imports
import argparse
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import math
import pandas as pd
import random

parser = argparse.ArgumentParser()
parser.add_argument("protein", help="Name of protein sequence", type=str)
parser.add_argument("up", help="Up: 30bp 5' of sequence", type=str)
parser.add_argument("wt", help="WT: sequence of interest", type=str)
parser.add_argument("down", help="Down: 30bp 3' of sequence ", type=str)
parser.add_argument("--vect", help="Vect: sequence of interest in vector", type=str, default=False)
parser.add_argument("--homo_len", help="Length of homology arm in fwd primer", type=int, default=20)
parser.add_argument("--melt_temp", help="Melting temp of fwd primer", type=int, default=50)
parser.add_argument("--rev_melt_temp", help="Melting temp of rev primer", type=int, default=55)
parser.add_argument("--syn_snp", help="Percentage of synonymous SNPs 0-1", type=float, default=.05)
parser.add_argument("--out_dir", help="Local output directory e.g. 'data/'", type=str, default=.05)

args = parser.parse_args()

# INPUT:
# sequence data
"""up = 'cgttggcacccagatttgaccttcctgacatgaaagaaacaaagtatact'.upper()
down = 'ggccaagttttcaaagcaaaacacagaattgacggaaagacttacgttat'.upper()

# FAKE WT PKR
wt = 'GTGGACAAGAGGTTTGGCATGGATTTTAAAGAAATAGAATTAATTGGCTCAGGTGGATTT'.upper()

# MCp74 window 1
vect = 'gtggacaagaggtttggcatggattttaaagaaatagaattaattggctcaggtggattt'.upper()

# PARAMETERS:
homo_len = 20
melt_temp = 50
rev_melt_temp = 55
protein = 'pkr'"""

# Assign args
protein = args.protein
up = args.up.upper()
wt = args.wt.upper()
down = args.down.upper()
if args.vect == False:
    vect = wt.upper()
else:
    vect = args.vect.upper()
homo_len = args.homo_len
melt_temp = args.melt_temp
rev_melt_temp = args.rev_melt_temp

# PREP:
# full wt and vector sequences
full_wt = up + wt + down
full_vect = up + vect + down

# full missense WT sequence
wt_codons = [wt[i:i+3] for i in range(0, len(wt), 3)]
mis_codons = wt_codons.copy()
mis_df = pd.read_csv('missense_table.csv')
mis_dict = pd.Series(mis_df['iupac'].values, index=mis_df['codon']).to_dict()
full_mis = up + ''.join([mis_dict[i] for i in wt_codons]) + down

# full vector synonymous sequence
vect_codons = [vect[i:i+3] for i in range(0, len(vect), 3)]
syn_codons = vect_codons.copy()
syn_df = pd.read_csv('syn_table.csv')
syn_dict = pd.Series(syn_df['iupac'].values, index=syn_df['codon']).to_dict()
full_syn = up + ''.join([syn_dict[i] for i in wt_codons]) + down

# START LOOP:
win_list = []
homo_list = []
rev_primer_list = []
fwd_primer_list = []
mut_windows_list = []
mut_name_list = []
full_window_list = []
full_name_list = []
full_primer_list = []

win = 1
index = len(up)
while index < len(up+wt):

    win_list.append(win)

    # 1. homology arm
    homo = full_vect[index-homo_len:index]
    homo_list.append(homo)

    # 2. reverse vector primer
    rev_temp = str(Seq(full_vect[:index]).reverse_complement())
    rev_primer = rev_temp[:15]
    while mt.Tm_NN(rev_primer) < rev_melt_temp:
        rev_primer = rev_temp[:len(rev_primer)+1]
    rev_primer_list.append(rev_primer)

    # 3. OPTIMIZED foward primer
    primer_end = index + (60-homo_len)
    if primer_end > len(full_vect):
        primer_end = len(full_vect)

    primer_start = primer_end - 15
    fwd_primer = full_vect[primer_start:primer_end]

    while mt.Tm_NN(fwd_primer) < melt_temp:
        primer_start -= 1
        fwd_primer = full_vect[primer_start:primer_end]

    # check if the primer is > 28bp
    if len(fwd_primer) > (60-homo_len-12):
        # fix mut window to 12, make a long primer
        primer_start = index + 12
        primer_end = primer_start + 15
        fwd_primer = full_vect[primer_start:primer_end]
        while True:
            fwd_primer = full_vect[primer_start:primer_end]
            if mt.Tm_NN(fwd_primer) > 50 and fwd_primer.upper().count('G') + fwd_primer.upper().count('C') > 8:
                break
            else:
                primer_end += 1

    else:
        # add or subtract a bp from the fwd primer to get mut_window in frame
        if ((primer_start)-index)%3 == 2:
            primer_start += 1
            fwd_primer = full_vect[primer_start:primer_end]

        elif ((primer_start)-index)%3 == 1:
            primer_start -= 1
            fwd_primer = full_vect[primer_start:primer_end]

    if primer_start > len(up+wt):
        # make the last fwd primer
        primer_start = len(up+wt)
        primer_end = primer_start+15
        fwd_primer = full_vect[primer_start:primer_end]
        while mt.Tm_NN(fwd_primer) < melt_temp:
            primer_end += 1
            fwd_primer = full_vect[primer_start:primer_end]
    fwd_primer_list.append(fwd_primer)

    # 4. mut window
    mut_len = (primer_start)-index
    mut_end = index+mut_len

    def codons_list(seq):
        return [seq[i:i+3] for i in range(0, len(seq), 3)]

    wt_list = codons_list(full_wt[index:mut_end])
    mis_list = codons_list(full_mis[index:mut_end])
    vect_list = codons_list(full_vect[index:mut_end])
    syn_win = full_syn[index:mut_end]

    mut_list = []
    mut_names = []
    for i,(wt_cod,mis_cod) in enumerate(zip(wt_list,mis_list)):
        for x in range(len(wt_cod)):
            if wt_cod[x] != mis_cod[x]: # make variant only if missense SNP accessible
                if random.random() < args.syn_snp:
                    mis_codon = wt_cod[:x]+"N"+wt_cod[x+1:]
                else:
                    mis_codon = wt_cod[:x]+mis_cod[x]+wt_cod[x+1:]
                mut_win = syn_win[:i*3] + mis_codon + syn_win[(i+1)*3:]
                mut_list.append(mut_win)
                mut_windows_list.append(mut_win)
                #mutation = f'{wt_cod[x]}-{mis_cod[x]}'
                mutation = f'{wt_cod}-{mis_codon}'
                mut_names.append(mutation)
            else:
                continue
    for i in mut_names:
        mut_name_list.append(i)

    # 5. Create full lists
    for i in mut_list:
        full_window_list.append(win)

    for i,mut in enumerate(mut_names):
        full_name_list.append(f'{protein}_fwd_{win}-{i+1}_{mut}')

    for i in mut_list:
        full_primer_list.append(homo+i+fwd_primer)

    index = primer_start
    win += 1

    if win > 20:
        break

len(mut_name_list)

rev_name_list = [f'{protein}_rev_{x}' for x in win_list]

# window summary df
df = pd.DataFrame(
    list(zip(win_list, fwd_primer_list, rev_primer_list)),
    columns=['window','fwd_primer','rev_primer'])
df['fwd_temp'] = df['fwd_primer'].apply(mt.Tm_NN).round(1)
df['fwd_gc'] = df['fwd_primer'].apply(GC).round(1)
df['rev_temp'] = df['rev_primer'].apply(mt.Tm_NN).round(1)
df['rev_gc'] = df['rev_primer'].apply(GC).round(1)

df.to_csv(f'{args.out_dir}{protein}_window-summary.csv', index=False)

df.head()

# doped primer summary df
df = pd.DataFrame(
    list(zip(full_window_list, full_name_list, full_primer_list)),
    columns=['window','name','full_oligo'])

df['full_oligo_len'] = df['full_oligo'].str.len()
df['homology_arm'] = df['window'].map(dict(zip(win_list,homo_list)))
df['homology_len'] = df['homology_arm'].str.len()
df['mut_window'] = df['name'].map(dict(zip(full_name_list,mut_windows_list)))
df['mut_len'] = df['mut_window'].str.len()
df['fwd_primer'] = df['window'].map(dict(zip(win_list,fwd_primer_list)))
df['fwd_temp'] = df['fwd_primer'].apply(mt.Tm_NN).round(1)
df['fwd_gc'] = df['fwd_primer'].apply(GC).round(1)
df['rev_primer'] = df['window'].map(dict(zip(win_list,rev_primer_list)))
df['rev_temp'] = df['rev_primer'].apply(mt.Tm_NN).round(1)
df['rev_gc'] = df['rev_primer'].apply(GC).round(1)

df.to_csv(f'{args.out_dir}{protein}_primer-summary.csv', index=False)

# idt file, sequence name, sequence
df_1 = pd.DataFrame(
    list(zip(full_name_list, full_primer_list)),
    columns=['sequence_name','sequence'])
df_2 = pd.DataFrame(
    list(zip(rev_name_list, rev_primer_list)),
    columns=['sequence_name','sequence'])
df = pd.concat([df_1, df_2])
df.reset_index(drop=True)

df.to_csv(f'{args.out_dir}{protein}_idt.csv', index=False)

# primer fasta file
f = open(f'{args.out_dir}{protein}_primers.fa', "w")
for x,y in zip(full_name_list,full_primer_list):
    f.write(f'>{x}\n')
    f.write(f'{y}\n')
for x,y in zip(rev_name_list,rev_primer_list):
    f.write(f'>{x}\n')
    f.write(f'{y}\n')
f.close()
