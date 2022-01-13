'''
Author: Anne de Jong 
date: June 2021
PyVersion: 3.7.6

Predict promoters in any DNA using best model
1) Best model selected on the basis of GC%

e.g.

sessiondir=/data/pg-bactprom/anne/00.PredictionOnly01
query=query
modelName=/data/pg-bactprom/anne/Parageobacillus_thermoglucosidasius.cnn_lstm_71/cnn_lstm.h5
python3 /data/pg-bactprom/scripts/ppp_Prediction_Only_Anne.py -sessiondir $sessiondir -query $query -modelName $modelName

'''
import os
import sys
import re
import argparse
import pandas as pd
import numpy as np
import D3GB ## Create a genome browser from fasta and gff file
from statistics import median_high
from tensorflow.keras.models import load_model
from scipy.cluster.hierarchy import ward, fcluster
from scipy.cluster.hierarchy import fclusterdata
from scipy.spatial.distance import pdist
from sklearn.metrics import confusion_matrix




# ---------------------------------------------------------------- parse parameters -------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Promoter Prediction')
parser.add_argument('-sessiondir', dest='sessiondir', help='Session Dir', nargs='?', default='.')
parser.add_argument('-query', dest='query', help='Full path and name genome basename;  without .fna .gff and .promoters.bed')
parser.add_argument('-modelName',dest='modelName', help='Full path and filename of the model', nargs='?', default='./cnn_lstm.h5')  # NEW
parser.add_argument('-promlen',dest='PromLen', help='Length of Promoter', nargs='?', default=71)
parser.add_argument('-prime',dest='prime', help='5 prime region after TSS', nargs='?', default=0)
parser.add_argument('-pval',dest='pvalue', help='p-value cutoff for initial prediction', nargs='?', default=0.99)
parser.add_argument('-out',dest='outPrefix', help='Prefix for Output files', nargs='?', default='ppp')
parser.add_argument('--version', action='version', version='Anne de Jong, version 2.0, Jan 2021')
args = parser.parse_args()


''' Uncomment for local use 
args.sessiondir = 'G:\\My Drive\\WERK\\PromoterPrediction\\Scripts_Anne\\results'
args.query = 'G:\My Drive\WERK\PromoterPrediction\Scripts_Anne\ppp_test'
args.modelName = "G:\My Drive\WERK\PromoterPrediction\Scripts_Anne\CNN_model.1000.h5"
args.modelName = "G:\My Drive\WERK\PromoterPrediction\Scripts_Anne\cnn_model.h5"
args.PromLen  = 51
args.prime  = 0
args.outPrefix = 'ppp'
args.pvalue = 0.99
'''

#  ---- be sure Arguments are handled as int or real value ------
args.PromLen  = int(args.PromLen)
args.prime  = int(args.prime)
args.pvalue = float(args.pvalue)

# ----- For clustering -----------------------------
MIN_CLUSTER_SIZE = 1
WINDOW_SIZE = 10
PROBABILITY_CUTOFF = 0.99

# ------ For analysis ------------------------------
SHIFT_FOR_NEXT_WINDOW = 1

#  gff header for genome annotation and  promoter files as BED
gff_header = ["chrom","db", "type", "start", "end", "name", "strand", "score", "description"]
bed_header = ["chrom","start", "end", 'name','score', "strand"]




''' ==============================  DEFINING FUNCTIONS =============================================== '''


def Anne_one_hot_encode(seq):
    mapping = dict(zip("ACGT", range(4)))
    seq2 = [mapping[i] for i in seq]
    return np.eye(4)[seq2]


def remPrevFormatting(fna_file):
    f = open(fna_file)
    output = ''
    for line in f:
        if line[0] != ">":
            output = output + line.strip()
    for i in range(0, len(output)):
        if output[i] == 'N' or output[i] == 'n':
                output = output.replace(output[i], 'G')
    return output

def getCleanSeq(fna_file):
	# get DNA and replace N or n by G or g to prevent error in training; G is considered as most save replacment
	DNA = ''
	with open(fna_file) as lines:
		for line in	lines: 
			if line[0] != ">": 
				DNA += line.strip()
	DNA = re.sub('[BDEFHIJKLMNOPQRSUVWXYZ]','G', DNA.upper())
	return DNA	

def getFastaKey(fna_file):
	# results the header between > and the first space
	fline=open(fna_file).readline().rstrip()
	keys = fline.split()
	return keys[0][1:]	



def makeQuerySet(DNA, window_shift):
	query_set = []
	for i in range(0, len(DNA)-args.PromLen-args.prime, window_shift):
		query_set.append(DNA[i:i+args.PromLen+args.prime])
	return query_set


def GCpercent(DNA):
	DNA = DNA.upper()
	return 100 * (DNA.count('G') + DNA.count('C')) / len(DNA)


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def write_log(S):
	# write to console and logfile
	print(S)
	f = open(args.sessiondir+'/00.ppp.log', "a")
	f.write(S + '\n')
	f.close()	






write_log('DNA file ='+args.query+'.fna')
write_log('GFF file ='+args.query+'.gff')
write_log('model file ='+args.modelName)
write_log('outfile = correct_predictions.txt')
write_log('================================================================')

	

''' ==============================  LOAD MODEL ========================================================== '''
write_log('Loading model and formatting DNA for prediction...')
model = load_model(args.modelName)


''' ==============================  LOAD DNA SEQUENCE ========================================================== '''

DNA_sense = getCleanSeq(args.query+'.fna')
DNA_antisense = reverse_complement(DNA_sense)
write_log('Number of bases in Clean DNA sequence: '+str(len(DNA_sense)))

test_sense_sequences     = makeQuerySet(DNA_sense,  1)
test_antisense_sequences = makeQuerySet(DNA_antisense,  1)

write_log("Encode sense strand")
input_sense_features = []
for sequence in test_sense_sequences:
	input_sense_features.append(Anne_one_hot_encode(sequence))


write_log("Encode anti-sense strand")
input_antisense_features = []
for sequence in test_antisense_sequences:
	input_antisense_features.append(Anne_one_hot_encode(sequence))


write_log('Stacking features')
np.set_printoptions(threshold=40)
sense_features     = np.stack(input_sense_features)
antisense_features = np.stack(input_antisense_features)

write_log('Sense and anti-sense-feature ready for prediction')



''' ==============================  MAKING PREDICTIONS ====================================================== '''
write_log(' ==============================  MAKING PREDICTIONS ====================================================== ')

write_log('Create pandas DataFrame as container for promoter data')
promoter_db = pd.DataFrame(columns=['position','score','strand','sequence'])
promoter_db.set_index('position')


# sense strand
write_log('Prediction on sense strand')
predicted_sense_labels = model.predict(sense_features)  # make a numpy.ndarray; original label and predicted label

# save the first 1000 for evaluation
#np.savetxt(args.sessiondir+'/'+"predicted_sense_labels.csv", predicted_sense_labels[0:1000], delimiter="\t")
write_log('\tGet promoters from sense strand with prediction pvalue > '+ str(args.pvalue))
predicted_sense_promoter_list = []
probabilityValueSense = []
for i in range(0,len(predicted_sense_labels)):
	if (predicted_sense_labels[i][1]) > args.pvalue:
		probabilityValueSense.append(str(predicted_sense_labels[i][1]))
		predicted_sense_promoter_list.append(test_sense_sequences[i])  # Get the DNA sequence
		promoter_db = promoter_db.append({'position':  i,'score':predicted_sense_labels[i][1], 'strand' :'+', 'sequence': test_sense_sequences[i]}, ignore_index=True )
write_log('\tNumber of promoters sense strand      : ' + str(len(predicted_sense_promoter_list)))


# anti-sense strand
write_log('Prediction on anti-sense strand')
predicted_antisense_labels = model.predict(antisense_features)  # make a numpy.ndarray; original label and predicted label

write_log('\tGet promoters from anti-sense strand with prediction pvalue > '+ str(args.pvalue))
predicted_antisense_promoter_list = []
probabilityValueAntisense = []
for i in range(0,len(predicted_antisense_labels)):
	if (predicted_antisense_labels[i][1]) > args.pvalue:
		probabilityValueAntisense.append(str(predicted_antisense_labels[i][1]))
		predicted_antisense_promoter_list.append(test_antisense_sequences[i])   # Get the DNA sequence
		promoter_db = promoter_db.append({'position': len(DNA_sense) - i,'score':predicted_antisense_labels[i][1], 'strand' :'-', 'sequence': test_antisense_sequences[i]}, ignore_index=True )
write_log('\tNumber of promoters anti-sense strand : ' + str(len(predicted_antisense_promoter_list)))

write_log('Total number of predicted promoters: ' + str(len(promoter_db)))
write_log('================================================================')



if len(promoter_db) == 0:
	write_log('No promoters found')
	write_log('================================================================')
	# write empty result file
	pd.DataFrame(columns=["No promoters found"]).to_csv(args.sessiondir+'/'+args.outPrefix+'.Promoters.txt')
	os._exit(1)		   

if len(promoter_db) >50000:
	write_log('Error; Too many promoters found. Exceeding limit of 50,000. Data size= '+str(len(promoter_db)))
	# write empty result file
	pd.DataFrame(columns=["Error; Too many promoters found"]).to_csv(args.sessiondir+'/'+args.outPrefix+'.Promoters.txt')
	os._exit(1)	




''' ==============================  PERFORM HIERARCHICAL CLUSTERING ========================================= '''


write_log(' ==============================  PERFORM HIERARCHICAL CLUSTERING ========================================= ')
# Clustering written by Daniel Kaptijn, modified by Anne
# For debugging: Optional reload prediction results DataFrame promoter_db
# For debugging: promoter_db_header=['position','score','strand','sequence']
# For debugging: promoter_db = pd.read_csv(args.sessiondir+'/'+args.outPrefix+'.Promoters.initial.txt',  sep ='\t' ,header=0,  names=promoter_db_header)


if len(predicted_sense_promoter_list) > 4:
	write_log('Cluster Sense Strand')
	Xs = promoter_db[['position','strand']].copy()	
	Xs = Xs[Xs['strand']=='+']
	Xs = Xs.drop('strand', 1)
	Xs['2D'] = 0 
	Zs = ward(pdist(Xs))
	sense_pred = fclusterdata(Xs, t=WINDOW_SIZE, criterion='distance')

	sense_dict = {}
	for i in range(0, len(sense_pred)):
		key = int(sense_pred[i])       # cluster ID
		value = int(Xs['position'][i])  # TSS position
		if key not in sense_dict.keys():
			sense_dict[key] = [value]
		else:
			new_value = [i for i in sense_dict[key]]
			new_value.append(value)
			sense_dict[key] = new_value

	for i in sense_dict.keys():
		if len(sense_dict[i]) >= MIN_CLUSTER_SIZE:
			cluster_centre = median_high(sense_dict[i])
			promoter_db.loc[(promoter_db['position'] == cluster_centre) & (promoter_db['strand'] == '+'), 'centre'] = cluster_centre


if len(predicted_antisense_promoter_list) > 4:
	write_log('Cluster Anti-Sense Strand')
	Xs = promoter_db[['position','strand']].copy()
	Xs = Xs[Xs['strand']=='-']
	Xs = Xs.drop('strand', 1)
	Xs = Xs.iloc[::-1].reset_index(drop=True)		# reverse order and reindex for anti-sense strand
	Xs['2D'] = 0
	Zs = ward(pdist(Xs))
	antisense_pred = fclusterdata(Xs, t=WINDOW_SIZE, criterion='distance')


	antisense_dict = {}
	for i in range(0, len(antisense_pred)):
		key = int(antisense_pred[i])
		value = int(Xs['position'][i])
		if key not in antisense_dict.keys():
			antisense_dict[key] = [value]
		else:
			new_value = [i for i in antisense_dict[key]]
			new_value.append(value)
			antisense_dict[key] = new_value

	for i in antisense_dict.keys():
		if len(antisense_dict[i]) >= MIN_CLUSTER_SIZE:
			cluster_centre = median_high(antisense_dict[i])
			promoter_db.loc[(promoter_db['position'] == cluster_centre) & (promoter_db['strand'] == '-'), 'centre'] = cluster_centre



write_log('\n')
write_log('Initial number of predictions:           ' +str(len(promoter_db)))
write_log('Number of predictions after clustering:  ' +str(promoter_db.centre.count()))
write_log('Saving clustered data: '+ args.outPrefix+'.Promoter.db.clustered.txt')

#promoter_db.to_csv(args.sessiondir+'/'+args.outPrefix+'.Promoter.db.clustered.txt', index = False, sep ='\t')


''' ==============================  ADD TSS POSITION ========================================= '''
# Anne, TSS added: The position is the start of the promoter. TSS is position + length of promoter
write_log('ADD TSS POSITION ')
promoter_db.loc[promoter_db['strand'] == '+', 'TSS'] = promoter_db['position'] + args.PromLen 
promoter_db.loc[promoter_db['strand'] == '-', 'TSS'] = promoter_db['position'] - args.PromLen 




''' ==============================  CLASSIFY PROMOTERS =================================================================== '''
# Anne, Classification added: Use the original features annotation (from gff) to map and classify TSS position
write_log('CLASSIFY PROMOTERS ')
genes_gff = pd.read_csv(args.query+".gff", header=None,  comment='#',sep='\t', names=gff_header)
convert_dict = { 'start': int, 'end': int }
genes_gff = genes_gff.astype(convert_dict)  # be sure that start and end are integers

# Here we only use the genes 
promoter_db['class'] = 'intergenic'
for index, feature in genes_gff.loc[genes_gff['type'] == 'gene'].iterrows():
	promoter_db.loc[(promoter_db['TSS'] > feature.start) & (promoter_db['TSS'] < feature.end), 'class'] = 'feature'


''' ==============================  EXPORT RESULT ======================================================================= '''


# Export the promoters
promoter_db['score'] = promoter_db['score'].apply(lambda x: round(x*10-9, 2))
promoter_db_centre = promoter_db.dropna()
del promoter_db_centre['centre']
promoter_db_centre.sort_values(by="TSS", ascending=True, inplace=True)
promoter_db_centre.reset_index(drop=True, inplace=True)


FastaKey = getFastaKey(args.query+'.fna')

# EXPORT promoters as GFF 
promoter_gff = pd.DataFrame(columns=gff_header)
df_row = pd.Series()
for index, row in promoter_db_centre.iterrows():
	df_row['chrom'] = FastaKey
	df_row['db'] = 'ProPr'
	df_row['type'] = 'promoter'
	df_row['strand'] = row['strand']
	if (row['strand'] == '+'):
		df_row['start'] = row['position']
		df_row['end']   = row['TSS']
	else:
		df_row['start'] = row['TSS']
		df_row['end']   = row['position']
	df_row['name'] = 'prom_'+format(index, '04d')
	df_row['score'] = row['score']
	df_row['description'] = 'ID='+df_row['name']+';TSS='+str(row['TSS'])+';class='+row['class']+';seq='+row['sequence']
	promoter_gff = promoter_gff.append(df_row, ignore_index=True)
promoter_gff.sort_values(by=['start']).to_csv(args.sessiondir+'/'+args.outPrefix+'.Promoters.gff', index = False, sep ='\t', columns=gff_header)

gffs = [genes_gff, promoter_gff]
pd.concat(gffs).sort_values(by=['start']).to_csv(args.sessiondir+'/'+args.outPrefix+'.Genes_Promoters.gff', index = False, sep ='\t', columns=gff_header,  header=False)

	

# EXPORT promoters as TABLE	
table_header= ["chrom","ID","TSS","strand","score","class","sequence","start","end"]
promoter_table = pd.DataFrame(columns=table_header)
df_row = pd.Series()
for index, row in promoter_db_centre.iterrows():
	df_row['chrom'] = FastaKey
	df_row['strand'] = row['strand']
	if (row['strand'] == '+'):
		df_row['start'] = row['position']
		df_row['end']   = row['TSS']
	else:
		df_row['start'] = row['TSS']
		df_row['end']   = row['position']
	df_row['ID'] = 'prom_'+format(index, '04d')
	df_row['TSS'] = row['TSS']
	df_row['score'] = row['score']
	df_row['sequence'] = row['sequence']
	df_row['class'] = row['class']
	promoter_table = promoter_table.append(df_row, ignore_index=True)
promoter_table.sort_values(by=['TSS']).to_csv(args.sessiondir+'/'+args.outPrefix+'.Promoters.txt', index = False, sep ='\t', columns=table_header)



''' ==============================  Make graphics ======================================================================= '''

print('D3GB folder = '+args.sessiondir+'/D3GB')
print('D3GB fasta  = '+args.query+'.fna')
print('D3GB gff    = '+args.sessiondir+'/'+args.outPrefix+'.Genes_Promoters.gff')
gb = D3GB.genomebrowser((args.query+'.fna',),directory = args.sessiondir+'/D3GB')
gb.addSequence(args.query+'.fna')
gb.addGFF(args.sessiondir+'/'+args.outPrefix+'.Genes_Promoters.gff')



