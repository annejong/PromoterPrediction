'''
Author: Anne de Jong 
date: Sept 2023

'''

import argparse
import pandas as pd


# ---------------------------------------------------------------- parse parameters -------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Add Regulon to UniProt')
parser.add_argument('-uniprot', dest='uniprot', help='Full path and name of uniprot database')
parser.add_argument('-regulon', dest='regulon', help='Full path and name of regulon database; .faa')
parser.add_argument('--version', action='version', version='Anne de Jong, version 1.0, Sept 2023')
args = parser.parse_args()





#  load gff  
gff_header = ["chrom","db", "type", "start", "end", "name", "strand", "score", "description"]
#gff = pd.read_csv(args.gff, header=None,  comment='#',sep='\t', names=gff_header)
#convert_dict = { 'start': int, 'end': int }
#gff = gff.astype(convert_dict)  # be sure that start and end are integers
#gff_genes  = gff.loc[gff['type'] == 'gene']
#gff_genes.sort_values(by=['start']).to_csv(args.gff, index = False, sep ='\t', columns=gff_header, header=None)



