# -*- coding: utf-8 -*-
"""
Created Feb 05 2020

@author: Anne
"""

# module load Python/3.5.1-foss-2016a

# e.g.
# python3 /data/ppp/frank.py -gff /data/g2d_mirror_genbank/Streptococcus_pneumoniae_D39V/ASM300349v1_genomic.gff -fna /data/g2d_mirror_genbank/Streptococcus_pneumoniae_D39V/ASM300349v1_genomic.g2d.fna -size 100000 -out frank

import sys
import pandas as pd
import argparse 
import re
import random

# ---------------------------------------------------------------- parse parameters -------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Frank maker')
parser.add_argument('-gff', dest='gff_file', help='GFF filename')
parser.add_argument('-fna', dest='fna_file', help='Complete genome FNA filename')
parser.add_argument('-size',dest='FragSize', help='size of fragments', nargs='?', default='100')
parser.add_argument('-frag',dest='FragNum', help='number of fragments', nargs='?', default='4')
parser.add_argument('-out', dest='outfile', help='Output filename prefix; prefix.fna and prefix.gff', nargs='?', default='./')
parser.add_argument('--version', action='version', version='Anne de Jong, version 1.0, Feb 2020')
args = parser.parse_args()
print('gff ='+args.gff_file)
print('fna ='+args.fna_file)
print('FragSize ='+args.FragSize)
print('outfile  ='+args.outfile)
print('===========================================================')

# ---------------------------------------------------------------- functions -------------------------------------------------------------------------------------------

def random_fragments_from_fna():
	fragments=[]
	for key in sorted(fna):
		l=len(fna[key])
		N=int(args.FragNum)
		S=int(args.FragSize)
		for i in range(N):
			start = 1 + round(i*(l/N),0)
			rndstart = random.randint(start,round(start+(l/N)-S,0))
			seq=fna[key][rndstart:rndstart+S]
			record = {'key': key, 'start': rndstart, 'end': rndstart+S, 'seq': seq}
			fragments.append(record)
	return fragments		

def load_fna():
	fasta = {}
	key=''
	with open(args.fna_file) as lines:
		for line in	lines:
			items = re.match("^>(.*)", line)
			if items:
				key= items.group(1)
				fasta[key]=''
			else:
				fasta[key]+=line
	return fasta			
		
				
def write_fna(header, seq):
	with open(args.outfile+'.fna', 'w') as f:	f.write('>'+header+'\n'+seq)	
	f.close()

# ---------------------------------------------------------------- main -------------------------------------------------------------------------------------------

# 1. Make the fragments
fna=load_fna()
fragments=random_fragments_from_fna()

# 2. Read the nine columns gff file using pandas
gff = pd.read_csv(args.gff_file,sep='\t',header=None, dtype=str, comment="#", names=["genome","db", "type", "start", "end", "dot1", "strand", "dot2", "description"])
convert_dict = {'start': int, 'end': int } 
gff = gff.astype(convert_dict)

# 3. Merge the fragments and the associated gff
ConcatFragment='';
GFFlist=[]
for fragment in fragments: 
	print('Fragment: ' + str(fragment['start'])+' - '+str(fragment['end']))
	select = gff[gff['start'].between(fragment['start']+100,fragment['end']-100)]
	select['start'] = select['start'] - fragment['start'] + len(ConcatFragment);
	select['end']   = select['end'] - fragment['start'] + len(ConcatFragment);
	GFFlist.append(select)
	ConcatFragment+=fragment['seq']

# 4. Write Sequence in fasta format and export the new GFF file 
write_fna('Frank',ConcatFragment)
pd.concat(GFFlist).to_csv(args.outfile+'.gff', sep='\t', index=False)

# 	
print('=======================')
print(ConcatFragment)
print(pd.concat(GFFlist))

  