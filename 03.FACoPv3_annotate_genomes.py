"""
Anne de Jong
2023 Feb
python3 /data/FACoPv2/03.FACoPv2_annotate_genomes.py -db refseq

python3 /data/FACoPv2/03.FACoPv2_annotate_genomes_HABROK.py -dbdir /scratch/public/genomes/prokaryotes -db refseq

COPY from HABROK
For single genome command line use -webserver true -dbdir

/data/FACoPv2_genomes/refseq/Escherichia_coli_O25bH4-ST131_EC958_EC958v1/Escherichia_coli_O25bH4-ST131_EC958_EC958v1_genomic

python3 /data/FACoPv2/03.FACoPv2_annotate_genomes.py -webserver true -dbdir /data/g2d_mirror_genbank/Escherichia_coli_O25b:H4-ST131_EC958
/data/g2d_mirror_genbank/Escherichia_coli_O25b:H4-ST131_EC958

"""

import sys
import pandas as pd
import argparse 
import subprocess
import os
import glob
import re


parser = argparse.ArgumentParser(description='FACoPv2 Annotation of Genome')
parser.add_argument('-db', dest='db', help='genbank or refseq', nargs='?', default='')
parser.add_argument('-dbdir', dest='dbdir', help='Database folder', nargs='?', default='/data/FACoPv2_genomes')
parser.add_argument('-fs', dest='fs', help='FileSearch prefix string', nargs='?', default='')
parser.add_argument('-webserver', dest='webserver', help='true or false', nargs='?', default='false')
parser.add_argument('--version', action='version', version='Anne de Jong, version 2.0, Feb 2023')
args = parser.parse_args()

programdir = os.path.dirname(os.path.realpath(__file__)) # get the root of the package
TransTermHP            = 'transterm'
TransTermHPexptermdat  = 'expterm.dat' # depending on TransTermHP transterm -p expterm.dat
FACoPv2_classification = '/data/FACoPv2/04.FACoPv2_genome_classification.sh'
codon_table_file       = '/data/FACoPv2/codon_table.txt'
gff_header = ["chrom","db", "type", "start", "end", "name", "strand", "score", "description"]
gff_descriptions = ['ID','Name','locus_tag','old_locus_tag']
min_intergenic_len = 15
max_intergenic_len = 600
max_operon_intergenic = 150


def webserver_log(text):
	f = open(args.dbdir+'/sessionprogress', "a")
	f.write(text+'<br>')
	f.close


def readMultiFasta(fna_file):
	# Note: string adding to dict is slow, here we use an intermediate string; seq
	fasta = {}
	key= ''
	seq=''
	with open(fna_file) as f:
		for line in f:
			if line.startswith(">"):
				if (key != ''):   # chk previous key
					fasta[key] = seq.strip().upper()
					seq = ''
				items = re.match(">(.*?)\s", line)
				if items: key=items.group(1)
			else: seq += line.rstrip()
		if (key != ''): fasta[key] = seq.strip().upper()	# add last record	
	return fasta	

# ======================== functions ============================================	

def unzip(FileSearch):
	print('Unzip files: '+FileSearch)
	files=glob.glob(FileSearch)
	for file in files:
		#print(file)
		cmd = "gunzip "+file
		s = subprocess.run([cmd], shell=True, universal_newlines=True,capture_output=True).stdout

def FACoP(genomeName):
	cmd = "/data/gsea_pro/FACoP/01.GSEApro_annotate_genomes.sh "+genomeName
	s = subprocess.run([cmd], shell=True, universal_newlines=True,capture_output=True).stdout

def GFF_add_decription_entries(GFF):
	for description in gff_descriptions: GFF[description]=""  # add entries to GFF
	for index, row in GFF.iterrows():
		items=row['description'].split(";")
		for item in items:
			entry=item.split("=")
			if entry and (entry[0] in gff_descriptions) : # add description entry if entry name exists
				GFF.loc[index,entry[0]] = entry[1]
	return GFF	

def reverse_complement(dna):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	dna = re.sub('[BDEFHIJKLMNOPQRSUVWXYZ]','G', dna)  # only keep GATC
	return ''.join([complement[base] for base in dna[::-1]])


def get_sequence(FASTA, gene):
	seq = FASTA[gene['chrom']][gene['start']-1:gene['end']]
	if gene['strand'] == '-': seq = reverse_complement(seq)
	return seq

def translate(seq):
	protein = ""
	for i in range(0, len(seq), 3) : 
		codon = seq[i:i+3]
		if codon in CODONS: protein += CODONS[codon]
	return protein

def load_codons(codon_table_file):
	codons={}
	with open(codon_table_file) as f:
		for line in f:
			items = line.rstrip().split("\t")
			if items[1]: codons[items[0]] = items[1]
	return codons			

def clean_seq(seq):
	unwanted = "'[]/+.!@#$;:!*%)(&^~="
	return ''.join( c for c in seq.strip() if c not in unwanted )

def make_coord_file(coordsfile, chrom, size):
	text = 'fakegene5\t1\t2\t'+chrom+'\n'
	text+= 'fakegene3\t'+str(size)+'\t'+str(size-1)+'\t'+chrom+'\n'
	f = open(coordsfile, "w")
	f.write(text)
	f.close
	
def make_ptt_file(pttfile, chrom, GFF):
	header = ['Location','Strand','Length','PID','Gene','Synonym','Code','COG','Product'] 
	ptt = pd.DataFrame()
	ptt['Location']=GFF.loc[GFF['chrom'] == chrom][["start", "end"]].astype(str).apply("..".join, axis=1)
	ptt['Strand']  =GFF.loc[GFF['chrom'] == chrom]['strand']
	ptt['Length']  =GFF.loc[GFF['chrom'] == chrom]['end']-genes['start']
	ptt['PID']     ='-'
	ptt['Gene']    =GFF.loc[GFF['chrom'] == chrom]['locus_tag']
	ptt['Synonym'] =GFF.loc[GFF['chrom'] == chrom]['locus_tag']
	ptt['Code']    ='-'
	ptt['COG']     ='-'
	ptt['Product'] ='-'
	pttfile = path+'/'+chrom+'.ptt'
	ptt.to_csv(pttfile, index = False, sep ='\t', columns = header)

	
	
# ======================== main ============================================	


FileSearch=args.dbdir+'/'+args.db+'/'+args.fs+'*/*.gff'   # add FileSearch Prefix, e.g., A for all strains starting with an A
if (args.webserver == 'true'):
	FileSearch=args.dbdir+'/*.gff';

CODONS = load_codons(codon_table_file)

# 1. Gunzip files if needed
if (args.webserver == 'false'):
	#unzip(args.dbdir+'/'+args.db+'/*/*.gz')
	unzip(args.dbdir+'/'+args.db+'/'+args.fs+'*/*.gz')


# 2. Annotate genomes
files=glob.glob(FileSearch)
fileCount=0
for file in files:
	if (args.webserver == 'true'): webserver_log('Annotating '+file)
	print(file)
	fileCount+=1
	if os.path.getsize(file) > 10000000: continue    # skip large gff files, which should nt be in but NCBI is not perfect
	path=os.path.dirname(file)
	genomeName=os.path.basename(file).replace('.gff','')
	print("\tCount: "+str(fileCount) + '/' + str(len(files)))
	print("\tpath="+path)
	print("\tname="+genomeName)

	# 1. Read genome sequence data
	print("\tload FASTA")
	FASTA = readMultiFasta(path+'/'+genomeName+'.fna')
	if (os.path.exists(path+'/'+genomeName+'.g2d.fna') is False):
		print("\twrite FNA to G2D format")
		f = open(path+'/'+genomeName+'.g2d.fna', "w")
		for chrom in FASTA:  # write the old ptt file for each fasta entry . This is needed for the old, but good, TranstermHP
			f.write('>'+chrom+'\n'+FASTA[chrom]+'\n')  # use clean headers
		f.close
		
	print("\tload GFF")
	GFF   = pd.read_csv(file, header=None,  comment='#',sep='\t', names=gff_header)
	convert_dict = { 'start': int, 'end': int }
	GFF = GFF.astype(convert_dict)  # be sure that start and end are integers
	GFF_add_decription_entries(GFF)
	
	# 2. Write GFF as Table if not exists
	filename=path+'/'+genomeName+'.g2d.table'
	if (os.path.exists(filename) is False):
		print("\tGFF to table")
		header = ["chrom","db", "type", "start", "end", "strand", 'ID','Name','locus_tag','old_locus_tag' ]
		GFF.sort_values(by=['start']).to_csv(filename, index = False, sep ='\t', columns = header)


	genes  = GFF.loc[GFF['type'] == 'gene']
	if (genes.shape[0] < 2): continue   # no genes are found in the GFF file, continue to the next file

	# 3. Write genes, intergenic and proteins in G2D FASTA format if not exists
	if (os.path.exists(path+'/'+genomeName+'.g2d.intergenic.ffn') is False):
		print("\tFNN, FAA and intergenic to G2D format")
		fnn = ''
		faa = ''
		intergenic = ''
		prevGene = pd.DataFrame()
		first = True
		for index, gene in genes.sort_values(by=['chrom','start']).iterrows():
			seq  = get_sequence(FASTA, gene)
			prot = translate(seq)
			fnn += '>'+gene['locus_tag']+'\n'+seq+'\n'
			faa += '>'+gene['locus_tag']+'\n'+prot+'\n'
			if first: first = False
			else:
				if gene['chrom'] == prevGene['chrom']: 
					intergenicSeq = FASTA[gene['chrom']][prevGene['end']:gene['start']]
					if len(intergenicSeq) > min_intergenic_len and len(intergenicSeq) < max_intergenic_len:
						if prevGene['strand'] == '-': intergenic += '>'+prevGene['locus_tag']+'\n'+reverse_complement(intergenicSeq)+'\n'
						if gene['strand']     == '+': intergenic += '>'+gene['locus_tag']+'\n'+intergenicSeq+'\n'
			prevGene=gene
		f = open(path+'/'+genomeName+'.g2d.fnn', "w")
		f.write(fnn)
		f.close
		f = open(path+'/'+genomeName+'.g2d.faa', "w")
		f.write(faa)
		f.close
		f = open(path+'/'+genomeName+'.g2d.intergenic.ffn', "w")
		f.write(intergenic)
		f.close
		
	
	# 4. Transcription Terminator prediction
	if (args.webserver == 'true'): webserver_log('Transcription Terminator prediction')
	gff_term_filename=path+'/'+genomeName+'.g2d.transterm'
	if (os.path.exists(gff_term_filename+'x') is False):
		print("\tTranscription Terminator prediction")
		fnafilename = path+'/'+genomeName+'.g2d.fna'
		gff_term = pd.DataFrame(columns=gff_header)
		count=0
		for chrom in FASTA:  # write the old ptt file for each fasta entry . This is needed for the old, but good, TranstermHP
			print('\t\t'+chrom)
			# Make .coords file
			coordsfile = path+'/'+chrom+'.coords'
			make_coord_file(coordsfile, chrom, len(FASTA[chrom]))   # make coords file needed for the transterm program (path+'/'+chrom+'.coords')
			# Make .ptt file 
			pttfile = path+'/'+chrom+'.ptt'
			make_ptt_file(pttfile, chrom, GFF)
			# Predict Terminators using the ptt file (optionally the .coords file can also be used)
			cmd = TransTermHP+ ' -p '+TransTermHPexptermdat+' '+fnafilename + ' '+pttfile
			s = subprocess.run([cmd], shell=True, universal_newlines=True,capture_output=True).stdout
			row = pd.Series(dtype='string')
			for line in s.split('\n'):
				items = re.match("\s+(TERM \d+)\s+(\d+)\s+-\s+(\d+)\s+(.)\s+(.)\s+(\d+)\s+(-?\d+.?\d*)\s+(-?\d+.?\d*)", line)
				if items:
					count+=1
					if (items.group(4) == '+'):
						start = int(items.group(2))
						end = int(items.group(3))
						seq=FASTA[chrom][start-1:end-1].strip()
					else:
						start = int(items.group(3))
						end = int(items.group(2))
						seq=reverse_complement(FASTA[chrom][start-1:end-1])
					if not pd.isna(start):
						row['chrom'] = chrom
						row['db'] = 'TransTermHP'
						row['type'] = 'Terminator'
						row['start'] = start
						row['end'] = end
						row['name'] = '.'
						row['strand'] = items.group(4)
						row['score'] = items.group(6)
						description  = 'ID=TERM_'+"{:04d}".format(count)
						description += ';locus_tag=TERM_'+"{:04d}".format(count)
						description += ';hairpin_score='+items.group(7)
						description += ';tail_score='+items.group(8)
						description += ';sequence='+seq
						row['description'] = description
						# OLD gff_term = gff_term.append(row, ignore_index=True)
						# using the loc indexer, append the series to the end of the df
						#gff_term.loc[len(gff_term)] =  row  # NEW
						gff_term = pd.concat([gff_term, row.to_frame().T], ignore_index=True)
		gff_term.start = gff_term.start.astype(int)
		gff_term.end = gff_term.end.astype(int)
		gff_term.sort_values(by=['chrom','start']).to_csv(gff_term_filename, index = False, sep ='\t', columns=gff_header, header=True)


	# 5. Operon prediction depends on Transcription Terminator prediction 
	if (args.webserver == 'true'): webserver_log('Operon prediction')
	operonfilename = path+'/'+genomeName+'.g2d.FACoP.OPERONS'
	if (os.path.exists(operonfilename+'x') is False):
		print("\tOPERON prediction")
		operons = pd.DataFrame()
		TransTerm = pd.read_csv(gff_term_filename, comment='#',sep='\t')
		TransTerm['center'] = (TransTerm.start+TransTerm.end) / 2
		operon=1
		row = pd.Series(dtype='string')
		first = True
		for index, gene in genes.sort_values(by=['chrom','start']).iterrows():
			if first: first = False
			else:
				if   gene['chrom'] != prevGene['chrom']:   operon+=1                      # other chrom				
				elif gene['strand'] != prevGene['strand']: operon+=1                   # gene in other strand
				elif gene['start']-prevGene['end'] > max_operon_intergenic: operon+=1 # large gap
				else:                                                                   # terminator present
					dfSize = TransTerm[TransTerm['center'].between(prevGene['end'], gene['start'])].size
					if (dfSize>0): operon+=1 
			row['locus_tag']   = gene['locus_tag']
			row['operonID']    = 'operon_'+"{:04d}".format(operon)
			row['description'] = gene['chrom']+';'+gene['description']
			# OLD operons = operons.append(row, ignore_index=True)
			#operons.loc[len(operons)] =  row # NEW
			operons = pd.concat([operons, row.to_frame().T], ignore_index=True)

			prevGene = gene	
		operons.sort_values(by=['locus_tag']).to_csv(operonfilename, index = False, sep ='\t', header=True)
		# make a copy of the FACoP operon file for backwards compatibility
		cmd = 'cp '+operonfilename + ' ' + path+'/'+genomeName+'.g2d.operons'
		print(cmd)
		s = subprocess.run([cmd], shell=True, universal_newlines=True,capture_output=True).stdout

		
	
	# 6. FACoPv2 genome classification  GO, COG etc
	if (args.webserver == 'true'): webserver_log('FACoPv2 genome classification  GO, COG, KEGG, IPR, eggNOG')
	G2Dfilename=path+'/'+genomeName+'.g2d'
	if (os.path.exists(G2Dfilename+'.FACoP.COG') is False):  # just check if one of the classification is already done
		print("\tCLASSIFICATION")
		cmd = FACoPv2_classification + ' ' + G2Dfilename
		print(cmd)
		s = subprocess.run([cmd], shell=True, universal_newlines=True,capture_output=True).stdout
	

# TODO
# Promoters
#





