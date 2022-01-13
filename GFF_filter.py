'''
Author: Anne de Jong 
date: Sept 2021
PyVersion: 3.7.6

GFF filter for D3GB genome viewer
only genes and RNA will be exported

e.g.

sessiondir=/data/pg-bactprom/anne/00.PredictionOnly01
query=query
modelName=/data/pg-bactprom/anne/Parageobacillus_thermoglucosidasius.cnn_lstm_71/cnn_lstm.h5
python3 /data/pg-bactprom/scripts/ppp_Prediction_Only_Anne.py -sessiondir $sessiondir -query $query -modelName $modelName

'''

import argparse
import pandas as pd


# ---------------------------------------------------------------- parse parameters -------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='GFF fiter')
parser.add_argument('-gff', dest='gff', help='Full path and name of gff file')
parser.add_argument('--version', action='version', version='Anne de Jong, version 1.0, Sept 2021')
args = parser.parse_args()


''' Uncomment for local use 
args.gff = 'C:\\Users\\Anne\\Downloads\\Staphylococcus_aureus_UA-S391_USA300_abscess_wound_isolate_ASM69587v1_genomic.gff'
'''

#  load gff  
gff_header = ["chrom","db", "type", "start", "end", "name", "strand", "score", "description"]
gff = pd.read_csv(args.gff, header=None,  comment='#',sep='\t', names=gff_header)
#convert_dict = { 'start': int, 'end': int }
#gff = gff.astype(convert_dict)  # be sure that start and end are integers


gff_genes  = gff.loc[gff['type'] == 'gene']

gff_genes.sort_values(by=['start']).to_csv(args.gff, index = False, sep ='\t', columns=gff_header, header=None)

