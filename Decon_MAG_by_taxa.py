#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
__author__ = 'Yaxin Xue'
__license__ = "GPL"
__version__ = "1.0.0"
__email__ = "yaxin.xue@uib.no, xue.ethan@gmail.com"

"""
import sys, os
import argparse
import pandas as pd
import Bio
from Bio import SeqIO

class MyParser(argparse.ArgumentParser):
   def error(self, message):
      sys.stderr.write('error: Check your Parameters!\n')
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

#%% functions
def init():
    example_text = '''Example:

    python Decon_MAG_by_taxa.py -r ./raw_MAG_fa/ -o ./cleaned_MAG_fa/ -l decon_mag_list.txt -m scaf2bin.txt -k kaiju_anno.txt -t 0.51

    '''
    parser = argparse.ArgumentParser(description='Extract dominant taxa of subset based on Kaiju annotation results\n',epilog=example_text, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser._optionals.title = 'Mandatory Arguments'
    parser.add_argument('-r', required=True, help='Path of your raw fasta folder')
    parser.add_argument('-o', required=True, help='Path of your output subset folder')
    parser.add_argument('-l', required=True, help='List of MAG ids you want to process')
    parser.add_argument('-m', required=True, help='Mapping information for Contig-MAG table')
    parser.add_argument('-k', required=True, help='Kaiju annotation result')
    parser.add_argument('-t', required=True, help='Threshold for dominant taxa percentage, default: 0.5', default=0.5, type=float)
    parser.print_help()
    return(parser)

def generate_ctg2mag_tax(kaiju, ctg2mag):
    # prase kaiju
    kaiju_df = pd.read_csv(kaiju, sep ='\t', header= None, skipinitialspace=True)
    taxa_info = kaiju_df.iloc[:, 3].str.split(';', expand=True)
    kaiju_df1 = pd.concat([kaiju_df.drop(3, axis=1), taxa_info],axis=1)
    kaiju_df1.columns = ['type','ctg_id','taxon_id','phylum','class','order', 'family', 'genus', 'species', 'other']
    # parse dastoolsca
    ctg2mag_pd = pd.read_csv(ctg2mag, sep = '\t', header = None, skipinitialspace=True)
    ctg2mag_pd.columns = ['ctg_id', 'bin_id']
    # merge info
    #kaiju_df1 = kaiju_df1.set_index('ctg_id',drop=False)
    #ctg2mag_pd = ctg2mag_pd.set_index('ctg_id', drop = False)
    ctg2mag_pd.index.name = None
    ctg2mag_pd_tax = ctg2mag_pd.merge(kaiju_df1, how='left', on='ctg_id')

    ctg2mag_wtax = ctg2mag_pd_tax
    return(ctg2mag_wtax)

def fill_NA_with_taxa(sel_ctg2mag_tax, taxa_rank):
    fillNA_ctg2mag = sel_ctg2mag_tax
    for idx, level in enumerate(taxa_rank):
        last_level = taxa_rank[idx - 1]
        if (level == 'phylum'):
            fillNA_ctg2mag.loc[:, level].fillna('P_na', inplace=True)
        else:
            fillNA_ctg2mag.loc[:, level].fillna(value = fillNA_ctg2mag.loc[:, last_level]+'_NA_'+level[0].upper(), inplace = True)
    return(fillNA_ctg2mag)

def calculate_taxa_perc(bin_id, taxa_rank, fillNA_ctg2mag):
    sel_ctg2mag_tax = fillNA_ctg2mag
    taxa_percent = dict()
    # ignore NA columns
    for level in taxa_rank:
        rank_count = sel_ctg2mag_tax.loc[:, level].value_counts(dropna = False)
        rank_perc = (rank_count/rank_count.sum()).sort_values(ascending = False)
        taxa_percent[level] = rank_perc
    return(taxa_percent)

def save_keeped_taxa(taxa_rank, taxa_percent, percent_cutoff):
    keeped_taxa = []
    for level in taxa_rank:
        keeped_rank = taxa_percent[level].loc[taxa_percent[level] > percent_cutoff]
        for taxa_name, taxa_perc in keeped_rank.items():
            keeped_name = '#'.join([level, str(taxa_name), str(round(taxa_perc,2))])
            keeped_taxa.append(keeped_name)
    return(keeped_taxa)


def print_cleaned_fa(raw_fa, cleaned_fa, keeped_ctg2mag):
    keeped_fa = []
    out_fa = open(cleaned_fa, 'w')
    for seq_record in SeqIO.parse(raw_fa, 'fasta'):
        seq_id = seq_record.id
        if seq_id in keeped_ctg2mag['ctg_id'].tolist():
            keeped_fa.append(seq_record)
    SeqIO.write(keeped_fa, out_fa, 'fasta')
    out_fa.close()
    return(keeped_fa)


def main():
    parser = init()
    args = parser.parse_args()
    # parse parameters
    raw_fa_path = args.r
    out_fa_path = args.o
    mag_list = args.l
    ctg2mag = args.m
    kaiju_anno = args.k
    tres_pct = float(args.t)
    if not os.path.exists(out_fa_path):
        os.mkdir(out_fa_path)
    #add taxa info to ctg
    ctg2mag_wtax = generate_ctg2mag_tax(kaiju_anno, ctg2mag)
    ctg2mag_wtax.to_csv('ctg2mag_wtax.csv', index=False)
    taxa_rank = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    for mag_id in open(mag_list):
        mag_id = mag_id.strip()
        print('Parse MAG: ', mag_id)
        sel_ctg2mag_wtax = ctg2mag_wtax.loc[ctg2mag_wtax['bin_id'] == mag_id]
        # replace space with underline for easy post process
        sel_ctg2mag_wtax = sel_ctg2mag_wtax.replace('\s+','_',regex=True)
        # step1: fill in NA with taxa
        fillNA_ctg2mag = fill_NA_with_taxa(sel_ctg2mag_wtax, taxa_rank)
        # step2: Calculate taxa percentage at each rank
        taxa_percent = calculate_taxa_perc(mag_id, taxa_rank, fillNA_ctg2mag)
        # step3: save output taxa above threshold
        output_taxa = save_keeped_taxa(taxa_rank, taxa_percent, tres_pct)
        # step4: print output subset fasta
        raw_fa = raw_fa_path+'/'+mag_id+'.fa'
        print(output_taxa)
        for output_taxa_info in output_taxa:
            (taxa_level, taxa_id, taxa_perc) = output_taxa_info.split('#')
            cleaned_name = mag_id +'.' + output_taxa_info.replace('#','__') + '.cleaned.fa'
            cleaned_fa = out_fa_path + '/' + cleaned_name
            output_ctg2mag = fillNA_ctg2mag[fillNA_ctg2mag[taxa_level] == taxa_id]
            output_fa = print_cleaned_fa(raw_fa, cleaned_fa, output_ctg2mag)
        print('Finish!\n')

if __name__ == "__main__":
   main()

#%%
