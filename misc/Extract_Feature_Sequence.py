# -*- coding: utf-8 -*-
"""
June 2021

@author: Anne

1) use linux to grep the genomes
grep  --include \*.gff -ie 'product=.*nifB' -r '/data/g2d_mirror/Azotobacter_vinelandii_DJ_DJ_ATCC_BAA-1303'
grep  --include \*.gff -ie 'product=.*bacteriocin' -r '/data/g2d_mirror/' > ExtractBacteriocins.sh
grep  --include \*.gff -ie 'product=.*lactose' -r '/data/g2d_mirror/Lactococcus' > Extractlactis.sh

2) use a downstream script (e.g. adapted from ExtractBacteriocins.sh) to process all grepped files: e.g. /data/ppp/ExtractNifB.sh
Run this script to extract upstream region of features
python3 /data/ppp/Extract_Feature_Sequence.py -genome /data/g2d_mirror/Acidianus_brierleyi_DSM_1651/ASM320183v2_genomic -feature NifB


"""
import sys
import re
import argparse
import pandas as pd


# ---------------------------------------------------------------- parse parameters -------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Promoter Prediction')
parser.add_argument('-sessiondir', dest='sessiondir', help='Session Dir', nargs='?', default='.')
parser.add_argument('-genome', dest='genome', help='Full path and name genome basename;  without .fna .gff and .promoters.bed')
parser.add_argument('-feature',dest='feature', help='Feature to extract')
parser.add_argument('-out',dest='outFile', help='Output file name', nargs='?', default='results.fna')
parser.add_argument('--version', action='version', version='Anne de Jong, version 2.0, Jan 2021')
args = parser.parse_args()


''' Uncomment for local use 
args.sessiondir = 'G:/My Drive/WERK/Python/Azotobacter'
args.genome = 'G:/My Drive/WERK/Python/Azotobacter/ASM38036v1_genomic'
args.feature = 'AVCA6_RS00030'
args.outFile = 'nifB.fna'
'''


def read_fasta(filename):
	results = {}
	dna=[]
	key=''
	with open(filename) as f:
		for line in f:
			if re.search("^>", line):
				if (len(dna)>0): results[key] = "".join(dna)
				my_list = re.match("^>(.*?)\s", line)
				key = my_list.group(1)
				results[key] = ''
			else:
				dna.append(line.rstrip())
	results[key] = "".join(dna)
	return results




def rev_compl(dna):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(dna))	


#  gff header for genome annotation and  promoter files as BED
gff_header = ["chrom","db", "type", "start", "end", "name", "strand", "score", "description"]
gff = pd.read_csv(args.genome+".gff", comment='#',header=None, sep='\t', names=gff_header)
gff = gff.astype({ 'start': int, 'end': int })  # be sure that start and end are integers

# load the genome sequence (can contain multiple entries)
fna = read_fasta(args.genome+".fna")

# get the feature
feature = gff.loc[(gff['type']=='CDS') & (gff['description'].str.contains(args.feature, case=True))]
print(args.genome, end='\t')
print(feature['chrom'].values[0], end='\t')
print(feature['start'].values[0], end='\t')
print(feature['end'].values[0], end='\t')
print(feature['strand'].values[0], end='\t')
print(feature['description'].values[0], end='\t')

# extract locus_tag from description
if (feature['strand'].values[0] == '+'):
	intergenic_start = feature['start'].values[0] - 150
	intergenic_end   = feature['start'].values[0] + 2
	DNA = fna[feature['chrom'].values[0]][intergenic_start:intergenic_end]

if (feature['strand'].values[0] == '-'):
	intergenic_start = feature['end'].values[0] - 3
	intergenic_end   = feature['end'].values[0] + 150
	DNA = rev_compl(fna[feature['chrom'].values[0]][intergenic_start:intergenic_end])
	
print(intergenic_start, intergenic_end, sep='\t', end='\t')
print(DNA,end='\t')
print('done')


