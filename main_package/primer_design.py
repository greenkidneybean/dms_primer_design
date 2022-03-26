#!/usr/bin/env python3
from . import codon_table
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

def homology_arm(seq_data, data_dict, args):
    start_index = data_dict['start_index']
    vector_seq = seq_data['vector_seq']

    homology_arm = vector_seq[start_index - args.homo_len:start_index] ### args.homo_len
    data_dict['homology_arm'] = homology_arm

    return data_dict

def reverse_primer(seq_data, data_dict, args):
    sub_window_name = data_dict['sub_window_name']
    start_index = data_dict['start_index']
    vector_seq = seq_data['vector_seq']

    reverse_seq = str(Seq(vector_seq[:start_index]).reverse_complement())
    reverse_primer = reverse_seq[:15]
    while mt.Tm_NN(reverse_primer) < args.rev_melt_temp: ### args.rev_melt_temp
        reverse_primer = reverse_seq[:len(reverse_primer)+1]
    data_dict['reverse_primer'] = reverse_primer

    reverse_primer_name = f'rev_{sub_window_name}'
    data_dict['reverse_primer_name'] = reverse_primer_name

    return data_dict

def forward_primer(seq_data, data_dict, args):
    start_index = data_dict['start_index']
    window_end = data_dict['window_end']
    vector_seq = seq_data['vector_seq']

    primer_end = start_index + (args.oligo_len - args.homo_len)
    if primer_end > window_end:
        primer_end == window_end

    primer_start = primer_end - 15
    forward_primer = vector_seq[primer_start:primer_end]

    while mt.Tm_NN(forward_primer) < args.melt_temp:
        primer_start -= 1
        forward_primer = vector_seq[primer_start:primer_end]

    # check if the primer is the max oligo length
    if len(forward_primer) > (args.oligo_len - args.homo_len - 12): # 12 is a minimum window size of 4 codons
        # fix mut window to 12, make a long primer
        primer_start = start_index + 12
        primer_end = primer_start + 15
        forward_primer = vector_seq[primer_start:primer_end]
        while True:
            forward_primer = vector_seq[primer_start:primer_end]
            if mt.Tm_NN(forward_primer) > args.melt_temp and forward_primer.upper().count('G') + forward_primer.upper().count('C') > 8:
                break
            else:
                primer_end += 1

    # even-out the primer length to accomodate codons
    else:
        # add or subtract a bp from the fwd primer to get mut_window in frame
        if (primer_start - start_index)%3 == 2:
            primer_start += 1
            forward_primer = vector_seq[primer_start:primer_end]

        elif (primer_start - start_index)%3 == 1:
            primer_start -= 1
            forward_primer = vector_seq[primer_start:primer_end]

    # making the last primer in a window
    if primer_start > window_end:
        primer_start = window_end
        primer_end = primer_start+15
        forward_primer = vector_seq[primer_start:primer_end]
        while mt.Tm_NN(forward_primer) < args.melt_temp:
            primer_end += 1
            forward_primer = vector_seq[primer_start:primer_end]

    data_dict['primer_start'] = primer_start
    data_dict['forward_primer'] = forward_primer

    return data_dict

def sub_window(seq_data, data_dict, args):
    primer_start = data_dict['primer_start']
    start_index = data_dict['start_index']
    window_end = data_dict['window_end']
    sub_window_name = data_dict['sub_window_name']
    wt_seq = seq_data['wt_seq']
    vector_seq = seq_data['vector_seq']
    gene_start = seq_data['gene_start']
    rng = seq_data['rng']

    # this may not work
    missense_dict, synonymous_dict, no_stop_dict, no_stop_syn_dict =  codon_table.iupac_codon_dicts()
    yeast_synonymous_dict = codon_table.synonymous_yeast_codons_dict()

    sub_window_len = (primer_start) - start_index
    sub_window_end = start_index + sub_window_len

    def codons_list(seq):
        return [seq[i:i+3] for i in range(0, len(seq), 3)]

    # removing mis_list and syn_list
    wt_list = codons_list(wt_seq[start_index:sub_window_end])
    vect_list = codons_list(vector_seq[start_index:sub_window_end])

    # generate synonymous vector codon list (top 2 codons for yeast)
    synonymous_win = [yeast_synonymous_dict[i].lower() for i in vect_list]

    # generate iupac missense codons list (with synonymous codons)
    iupac_codons = []
    syn_bool_list = []
    no_stop_list = []
    for i, wt_codon in enumerate(wt_list):
        syn_bool = rng.choice([True, False], p=[args.syn_snp_rate, 1-args.syn_snp_rate]) ### args.syn_snp_rate
        syn_bool_list.append(syn_bool)

        no_stop_bool = rng.choice([True, False], p=[args.stop_rate, 1-args.stop_rate]) ### args.stop_rate
        no_stop_list.append(no_stop_bool)

        # missense_dict, synonymous_dict, no_stop_dict, no_stop_syn_dict

        # need a sub-step to check if codon contains stop codon
        # need a function to check if the wt codon encodes a missense stop codon
        if syn_bool and no_stop_bool:
            # use no_stop_syn_dictionary
            iupac_codons.append(no_stop_syn_dict[wt_codon])
        elif syn_bool and not no_stop_bool:
            # use synonymous_dictionary
            iupac_codons.append(synonymous_dict[wt_codon])
        elif no_stop_bool and not syn_bool:
            # use no_stop_dict
            iupac_codons.append(no_stop_dict[wt_codon])
        else:
            # use missense dict
            iupac_codons.append(missense_dict[wt_codon])

    # make full-length oligo (homology arm, sub-window, primer), generate dataframe
    for i, iupac_list in enumerate(iupac_codons):
        aa_position = int((((start_index-gene_start)/3)+1)+i)
        # could enumerate this out to get the aas
        for iupac_codon in iupac_list:
            # get AAs encoded by iupac codon
            iupac_aa = codon_table.iupac_to_aa(iupac_codon)

            # place iupac_codon into sub_window
            sub_window = ''.join(synonymous_win[:i] + [iupac_codon] + synonymous_win[i+1:])

            codon_sub = wt_list[i] + str(aa_position) + iupac_codon
            forward_primer_name = f'{sub_window_name}_{codon_sub}'
            full_forward_primer = data_dict['homology_arm'] + sub_window + data_dict['forward_primer']

            # add values to data_dict
            dict_keys = ['name','codon_sub','wt','position','iupac', 'iupac_aa','sub_window', 'primer', 'synonymous_codons', 'no_stop_codons']
            dict_values = [forward_primer_name, codon_sub, wt_list[i], aa_position, iupac_codon, iupac_aa, sub_window, full_forward_primer, syn_bool_list[i], no_stop_list[i]]
            for (key,value) in zip(dict_keys,dict_values):
                data_dict[key] = value

            # append data_dict to dataframe
            seq_data['df'] = seq_data['df'].append(data_dict, ignore_index=True)

            # write primers to .fasta file
            seq_data['fasta_file'] = seq_data['fasta_file'] + [f">{forward_primer_name}\n", f"{full_forward_primer}\n"]

    return seq_data, data_dict
