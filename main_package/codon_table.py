#!/usr/bin/env python3
import pkg_resources
from itertools import product
from Bio.Seq import Seq
import pandas as pd
import itertools

iupac_dict = {'A':'A','C':'C','G':'G','T':'T','AC':'M','AG':'R','AT':'W','CG':'S','CT':'Y','GT':'K','ACG':'V','ACT':'H','AGT':'D','CGT':'B','ACGT':'N'}
rev_iupac_dict = {value:key for key,value in iupac_dict.items()}

def unique_missense_variants(codon, codon_table='Standard'):
    """Return all unique AA variants for a given codon
    INPUT: codon (str)
    RETURN: unique AAs (str)
    """
    aa = str(Seq(codon).translate(table=codon_table))
    nucs = 'ACTG'
    missense_aa_list = []
    for i in range(len(codon)):
        for n in nucs:
            new_codon = codon[:i] + n + codon[i+1:]
            new_aa = str(Seq(new_codon).translate(table=codon_table))
            if new_aa == aa:
                continue
            else:
                missense_aa_list.append(new_aa)
    return ''.join(set(missense_aa_list))

def contains_stop_missense_variant(codon, codon_table='Standard'):
    """Check if codon contains a stop missense variant
    INPUT: codon (str)
    RETURN: boolean (bool)
    """
    missense_aa = unique_missense_variants(codon, codon_table=codon_table)
    return "*" in missense_aa

def iupac_to_aa(iupac_codon):
    """Return string of AAs encoded by input iupac missense codon"""
    nuc_lists = [list(rev_iupac_dict[n]) for n in iupac_codon]
    codon_list = [''.join(i) for i in list(itertools.product(*nuc_lists))]
    aa_list = [str(Seq(codon).translate()) for codon in codon_list]
    return ''.join(aa_list)

def iupac_missense_codon_df(codon_table='Standard'):
    nucleotides = 'ACGT'
    iupac_dict = {'A':'A','C':'C','G':'G','T':'T','AC':'M','AG':'R','AT':'W','CG':'S','CT':'Y','GT':'K','ACG':'V','ACT':'H','AGT':'D','CGT':'B','ACGT':'N'}

    codon_list = []
    aa_list = []
    position_list = []
    nucleotides_list = []
    missense_codons_list = []
    missense_aa_list = []
    iupac_list = []
    iupac_codon_list = []

    for codon in list(''.join(w) for w in product(nucleotides, repeat=3)):
        aa = Seq(codon).translate(table=codon_table)[0]

        # loop through each position in codon
        for position in range(3):
            new_aas = []
            iupac_n = []
            new_codons = []
            # looping through each nucleotide
            for i in nucleotides:
                new_codon = codon[:position] + i + codon[position + 1:]
                new_aa = Seq(new_codon).translate(table=codon_table)[0]
                if not new_aa == aa:
                    new_aas.append(new_aa)
                    iupac_n.append(i)
                    new_codons.append(new_codon)
                else:
                    continue

            #iupac
            if not iupac_n: # check if iupac_n is empty
                iupac = codon[position]
            else:
                for i in iupac_dict.keys():
                    if set(i) == set(iupac_n):
                        iupac = iupac_dict[i]

            # make assignments here
            codon_list.append(codon)
            aa_list.append(aa)
            position_list.append(position)
            nucleotides_list.append(''.join(iupac_n))
            missense_codons_list.append(' '.join(new_codons))
            missense_aa_list.append(''.join(set(new_aas)))


            for i in iupac_dict.keys():
                if set(i) == set(iupac_n):
                    iupac = iupac_dict[i]
            iupac_list.append(iupac)
            iupac_codon_list.append(codon[:position] + iupac + codon[position + 1:])

    # create dictionary
    codon_level_dict = {
        'codon':codon_list,
        'aa':aa_list,
        'position':position_list,
        'missense_nucleotides':nucleotides_list,
        'missense_codons':missense_codons_list,
        'missense_aa':missense_aa_list,
        'iupac':iupac_list,
        'iupac_codon':iupac_codon_list
    }

    temp_df = pd.DataFrame.from_dict(codon_level_dict)
    return temp_df

def iupac_missense_codon_dict(codon_table='Standard'):
    temp_df = iupac_missense_codon_df(codon_table=codon_table)
    temp_dict = temp_df[temp_df.missense_aa != ''].groupby('codon')['iupac_codon'].apply(list).to_dict()
    return temp_dict

def iupac_synonymous_codon_df(codon_table='Standard'):
    nucleotides = 'ACGT'
    iupac_dict = {'A':'A','C':'C','G':'G','T':'T','AC':'M','AG':'R','AT':'W','CG':'S','CT':'Y','GT':'K','ACG':'V','ACT':'H','AGT':'D','CGT':'B','ACGT':'N'}

    codon_list = []
    aa_list = []
    position_list = []
    nucleotides_list = []
    synonymous_codons_list = []
    synonymous_aa_list = []
    iupac_list = []
    iupac_codon_list = []

    for codon in list(''.join(w) for w in product(nucleotides, repeat=3)):
        aa = Seq(codon).translate(table='Standard')[0]

        # loop through each position in codon
        for position in range(3):
            new_aas = []
            iupac_n = []
            new_codons = []
            # looping through each nucleotide
            for i in nucleotides:
                new_codon = codon[:position] + i + codon[position + 1:]
                new_aa = Seq(new_codon).translate(table='Standard')[0]
                if new_aa == aa:
                    new_aas.append(new_aa)
                    iupac_n.append(i)
                    new_codons.append(new_codon)
                else:
                    continue

            #iupac
            if not iupac_n: # check if iupac_n is empty
                iupac = codon[position]
            else:
                for i in iupac_dict.keys():
                    if set(i) == set(iupac_n):
                        iupac = iupac_dict[i]

            # make assignments here
            codon_list.append(codon)
            aa_list.append(aa)
            position_list.append(position)
            nucleotides_list.append(''.join(iupac_n))
            synonymous_codons_list.append(' '.join(new_codons))
            synonymous_aa_list.append(''.join(set(new_aas)))


            for i in iupac_dict.keys():
                if set(i) == set(iupac_n):
                    iupac = iupac_dict[i]
            iupac_list.append(iupac)
            iupac_codon_list.append(codon[:position] + iupac + codon[position + 1:])

    # create dictionary
    codon_level_dict = {
        'codon':codon_list,
        'aa':aa_list,
        'position':position_list,
        'synonymous_nucleotides':nucleotides_list,
        'synonymous_codons':synonymous_codons_list,
        'synonymous_aa':synonymous_aa_list,
        'iupac':iupac_list,
        'iupac_codon':iupac_codon_list
    }

    temp_df = pd.DataFrame.from_dict(codon_level_dict)
    return temp_df

def iupac_synonymous_codon_dict(codon_table='Standard'):
    # standard table will not have ATG and TGG synonymous iupac variants
    nucleotides = 'ACGT'
    temp_df = iupac_synonymous_codon_df(codon_table=codon_table)
    temp_df = temp_df[~temp_df.iupac.str.contains('|'.join(list(nucleotides)))]
    temp_dict = temp_df.groupby('codon')['iupac_codon'].apply(list).to_dict()
    return temp_dict


#### I need something that gets the absolute path of the .csv file to be imported

# SCRATCH
def selected_iupac_codons_dict():
    """return codon table of selected missense codons
    PROBLEM: codon table must be in directory where function is being called"""
    df = pd.read_csv('dms_codon_table_v2.csv')
    df.fillna('', inplace=True)
    sele_dict = df.query('sele_iupac_codon != ""').groupby('codon')['sele_iupac_codon'].apply(list).to_dict()
    for key,value in sele_dict.items():
        sele_dict[key] = list(itertools.chain.from_iterable([codon.split(' ') for codon in value]))
    return sele_dict

# SCRATCH
def synonymous_iupac_codons_dict():
    """return codon table including synonymous codons
    PROBLEM: codon table must be in directory where function is being called"""
    df = pd.read_csv('dms_codon_table_v2.csv')
    df.fillna('', inplace=True)
    syn_dict = df.query('syn_iupac_codon != ""').groupby('codon')['syn_iupac_codon'].apply(list).to_dict()
    for key,value in syn_dict.items():
        syn_dict[key] = list(itertools.chain.from_iterable([codon.split(' ') for codon in value]))
    return syn_dict


def iupac_codon_dicts():
    """Returns four mapping dictionaries to generate missense variants
    RETURNS:
    - missense_dict
    - synonymous_dict
    - no_stop_dict
    - no_stop_syn_dict
    """
    stream = pkg_resources.resource_stream(__name__, 'data/final_codon_table.csv')
    df = pd.read_csv(stream)
    df.fillna('', inplace=True)
    col_list = ['sele_iupac_codon', 'syn_iupac_codon', 'no_stop_iupac_codon', 'no_stop_syn_iupac_codon']

    # check that column names exist
    for col in col_list:
        if col not in df.columns.tolist():
            print(f"ERROR: Column '{col}' not contained in file '{codon_table}'")
            return
    dict_list = []
    for col in col_list:
        temp_dict = df.query(f'{col} != ""').groupby('codon')[col].apply(list).to_dict()
        for key,value in temp_dict.items():
            temp_dict[key] = list(itertools.chain.from_iterable([codon.split(' ') for codon in value]))
        dict_list.append(temp_dict)
    return dict_list

def synonymous_yeast_codons_dict():
    """Return mapping dictionary of doped synonymous codons optimized for yeast"""
    stream = pkg_resources.resource_stream(__name__, 'data/yeast_synonymous_codon_table.csv')
    df = pd.read_csv(stream)
    syn_dict = dict(zip(df.codon, df.iupac))
    return syn_dict
