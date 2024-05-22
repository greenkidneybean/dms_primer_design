#!/usr/bin/env python3

# functions to design DMS primers
# many functions require "args" input from script
import pkg_resources
import re
import yaml
from . import codon_table
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

def merge_config_files(custom_config, default_config=None):
    """Update two config dictionaries recursively.
    INPUT:
    custom_config: path to custom config file
    default_config: path to custom default file, optional
    """
    with open(custom_config) as cf_file:
        dict1 = yaml.safe_load(cf_file.read())
    if not default_config:
        default_config_file = pkg_resources.resource_filename(__name__, 'data/default_config.yaml')
        with open(default_config_file) as cf_file:
            dict2 = yaml.safe_load(cf_file.read())
    else:
        with open(default_config) as cf_file:
            dict2 = yaml.safe_load(cf_file.read())
    for k, v in dict2.items():
        if k not in dict1:
            dict1[k] = v
    return dict1

def homology_arm(seq_data, data_dict, cfg):
    start_index = data_dict['start_index']
    vector_seq = seq_data['vector_seq']

    homology_arm = vector_seq[start_index - cfg['homology_length']:start_index]
    data_dict['homology_arm'] = homology_arm

    return data_dict

def reverse_primer(seq_data, data_dict, cfg):
    sub_window_name = data_dict['sub_window_name']
    start_index = data_dict['start_index']
    vector_seq = seq_data['vector_seq']

    reverse_seq = str(Seq(vector_seq[:start_index]).reverse_complement())
    reverse_primer = reverse_seq[:15]
    while mt.Tm_NN(reverse_primer) < cfg['rev_primer_melt_temp']:
        reverse_primer = reverse_seq[:len(reverse_primer)+1]
    data_dict['reverse_primer'] = reverse_primer

    reverse_primer_name = f'rev_{sub_window_name}'
    data_dict['reverse_primer_name'] = reverse_primer_name

    # write reverse sub-window primer to .fasta file
    seq_data['fasta_file'] = seq_data['fasta_file'] + [f">{reverse_primer_name}\n", f"{reverse_primer}\n"]

    return data_dict

def forward_primer(seq_data, data_dict, cfg):
    start_index = data_dict['start_index']
    window_end = data_dict['window_end']
    vector_seq = seq_data['vector_seq']

    primer_end = start_index + (cfg['max_oligo_length'] - cfg['homology_length'])
    if primer_end > window_end:
        primer_end == window_end

    primer_start = primer_end - 15
    forward_primer = vector_seq[primer_start:primer_end]

    while mt.Tm_NN(forward_primer) < cfg['fwd_primer_melt_temp']:# and len(re.findall("[GC]", forward_primer.upper())) < 8: # check for GC count?
        primer_start -= 1
        forward_primer = vector_seq[primer_start:primer_end]

    # check if the primer is the max oligo length
    if len(forward_primer) > (cfg['max_oligo_length'] - cfg['homology_length'] - 12): # 12 is a minimum window size of 4 codons
        # fix mut window to 12, make a long primer
        primer_start = start_index + 12
        primer_end = primer_start + 15
        forward_primer = vector_seq[primer_start:primer_end]
        while True:
            forward_primer = vector_seq[primer_start:primer_end]
            if mt.Tm_NN(forward_primer) > cfg['fwd_primer_melt_temp'] and len(re.findall("[GC]", forward_primer.upper())) > 8:
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
    if primer_start > window_end or window_end - primer_start <= 3:
        primer_start = window_end
        primer_end = primer_start+15
        forward_primer = vector_seq[primer_start:primer_end]
        while mt.Tm_NN(forward_primer) < cfg['fwd_primer_melt_temp']:# and len(re.findall("[GC]", forward_primer.upper())) < 8:
            primer_end += 1
            forward_primer = vector_seq[primer_start:primer_end]

    # reduce primer if whole window is 61 bases
    if (primer_end - start_index + cfg['homology_length']) == 61:
        forward_primer = forward_primer[:-1]

    data_dict['primer_start'] = primer_start
    data_dict['forward_primer'] = forward_primer

    return data_dict

def sub_window(seq_data, data_dict, cfg):
    primer_start = data_dict['primer_start']
    start_index = data_dict['start_index']
    window_end = data_dict['window_end']
    sub_window_name = data_dict['sub_window_name']
    wt_seq = seq_data['wt_seq']
    vector_seq = seq_data['vector_seq']
    gene_start = seq_data['gene_start']
    rng = seq_data['rng']

    syn_var_rate = cfg['synonymous_variant_rate']
    remove_stop_rate = cfg['remove_stop_variant_rate']
    table = cfg['codon_table']

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

    # generate list of iupac missense codons to use
    # check to add synonymous variants and remove stop codons
    iupac_codons = []
    add_synonymous_codon_list = []
    contains_stop_list = []
    remove_stop_list = []
    for wt_codon in wt_list:
        # include synonymous variants (bool)
        syn_bool = rng.choice([True, False], p=[syn_var_rate, 1-syn_var_rate])
        add_synonymous_codon_list.append(syn_bool)

        # check if codon contains stop missense variants
        stop_bool = codon_table.contains_stop_missense_variant(wt_codon, table)
        contains_stop_list.append(stop_bool)

        # if codon contains stop variants
        if stop_bool:
            # remove stop variant (bool)
            remove_stop_bool = rng.choice([True, False], p=[remove_stop_rate, 1-remove_stop_rate])
        else:
            remove_stop_bool = False
        remove_stop_list.append(remove_stop_bool)

        # assign iupac codons for wt_codon
        if syn_bool and remove_stop_bool:
            # use no_stop_syn_dictionary, add syn and remove stops
            iupac_codons.append(no_stop_syn_dict[wt_codon])
        elif syn_bool and not remove_stop_bool:
            # use syn_dict, add syn and keep stops
            iupac_codons.append(synonymous_dict[wt_codon])
        elif not syn_bool and remove_stop_bool:
            # use the no_stop_dict, no syn and remove stops
            iupac_codons.append(no_stop_dict[wt_codon])
        else:
            # use missense_dict, no syn and keep stops
            iupac_codons.append(missense_dict[wt_codon])

    # make full-length oligo (homology arm, sub-window, primer), generate dataframe
    for i, iupac_list in enumerate(iupac_codons):
        aa_position = int((((start_index-gene_start)/3)+1)+i)
        # could enumerate this out to get the aas
        for iupac_codon in iupac_list:
            # get AAs encoded by iupac codon
            iupac_aa = codon_table.iupac_to_aa(iupac_codon)
            codon_subs = codon_table.codon_substitutions(wt_list[i], aa_position, iupac_codon)

            # place iupac_codon into sub_window
            sub_window = ''.join(synonymous_win[:i] + [iupac_codon] + synonymous_win[i+1:])

            iupac_sub = wt_list[i] + str(aa_position) + iupac_codon
            forward_primer_name = f'{sub_window_name}_{iupac_sub}'
            full_forward_primer = data_dict['homology_arm'] + sub_window + data_dict['forward_primer']

            # add values to data_dict
            dict_keys = [
                'name',
                'iupac_sub',
                'wt_codon',
                'position',
                'iupac_codon',
                'iupac_aa',
                'sub_window',
                'primer',
                'add_synonymous_codon',
                'contains_missense_stop',
                'remove_missense_stop_codon',
                'codon_subs'
                ]
            dict_values = [
                forward_primer_name,
                iupac_sub,
                wt_list[i],
                aa_position,
                iupac_codon,
                iupac_aa,
                sub_window,
                full_forward_primer,
                add_synonymous_codon_list[i],
                contains_stop_list[i],
                remove_stop_list[i],
                codon_subs
                ]
            for (key,value) in zip(dict_keys,dict_values):
                data_dict[key] = value

            # append data_dict to dataframe
            seq_data['df'] = seq_data['df'].append(data_dict, ignore_index=True)

            # write forward variant primer to .fasta file
            seq_data['fasta_file'] = seq_data['fasta_file'] + [f">{forward_primer_name}\n", f"{full_forward_primer}\n"]

    return seq_data, data_dict
