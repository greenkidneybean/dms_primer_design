#!/usr/bin/env python3

from itertools import product
from Bio.Seq import Seq
import pandas as pd
import itertools

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

def selected_missense_codons_dict():
    """return codon table of selected missense codons
    PROBLEM: codon table must be in directory where function is being called"""
    df = pd.read_csv('dms_codon_table_v2.csv')
    df.fillna('', inplace=True)
    sele_dict = df.query('sele_missense_codons != ""').groupby('codon')['sele_missense_codons'].apply(list).to_dict()
    for key,value in sele_dict.items():
        sele_dict[key] = list(itertools.chain.from_iterable([codon.split(' ') for codon in value]))
    return sele_dict

def synonymous_missense_codons_dict():
    """return codon table including synonymous codons
    PROBLEM: codon table must be in directory where function is being called"""
    df = pd.read_csv('dms_codon_table_v2.csv')
    df.fillna('', inplace=True)
    syn_dict = df.query('syn_missense_codons != ""').groupby('codon')['syn_missense_codons'].apply(list).to_dict()
    for key,value in syn_dict.items():
        syn_dict[key] = list(itertools.chain.from_iterable([codon.split(' ') for codon in value]))
    return syn_dict

def synonymous_yeast_codons_dict():
    """return iupac codon dictionary for synonymous codons most frequently used in yeast
    PROBLEM: codon table must be in directory where function is being called"""
    df = pd.read_csv('yeast_synonymous_codon_table.csv')
    temp_df = dict(zip(df.codon, df.iupac))
    return temp_df