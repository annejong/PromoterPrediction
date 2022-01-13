#!/bin/bash

# Main script for promoter prediction in prokaryotes
#
# Anne de Jong, June 2021 
#

PROGRAMDIR=/data/ppp

SESSIONDIR=$1
#SESSIONDIR='/data/ppp/misc'
ModelSelected=$2

# check code injection
COUNT1=`grep -rnw $SESSIONDIR -e '\?php' | wc -l `
COUNT2=`grep -rnw $SESSIONDIR -e '\?PHP' | wc -l `
if ( (($COUNT1+$COUNT2)) > 0)
then
	echo 'PHP code injection detected'
	echo 'PHP code injection detected' > $SESSIONDIR/00.ppp.log
	exit 1
else
	echo 'clean files'
fi	
 

# Get filenames
DNAFILES=`find $SESSIONDIR -type f \( -iname \*.fna -o -iname \*.fasta \)`
GFFFILES=`find $SESSIONDIR -type f \( -iname \*.gff -o -iname \*.gff3  \)`
COUNTDNA=`echo $DNAFILES | wc -w`

# Check if DNA files are found
COUNTDNA=`echo $DNAFILES | wc -w`
if [[ $COUNTDNA >0 ]]
then
	echo 'Number of DNA files found: ' $COUNTDNA >> $SESSIONDIR/00.ppp.log
else
	echo 'No valid DNA file found' >> $SESSIONDIR/00.ppp.log
	exit 1
fi	

# List the files
echo 'Files found in session folder;' >> $SESSIONDIR/00.ppp.log
for filename in $FILES;
do       
  echo $filename >> $SESSIONDIR/00.ppp.log
done;


# Combine all DNAFILES
mkdir $SESSIONDIR/query_combined_files
cat $DNAFILES > $SESSIONDIR/query_combined_files/query.fna

# Combine all GFFFILES
cat $GFFFILES > $SESSIONDIR/query_combined_files/query.gff

# For the D3GB viewer only use the type=gene lines to prevent errors in D3GB
python3 $PROGRAMDIR/GFF_filter.py -gff $SESSIONDIR/query_combined_files/query.gff

# Select model on the basis of GC content
case $ModelSelected in
	GC20|GC21|GC22|GC23|GC24|GC25|GC26|GC27|GC28|GC29) ModelSelected="FGC3040_30";;
	GC30|GC31) ModelSelected="FGC3040_30";;
	GC32|GC33) ModelSelected="FGC3040_32";;
	GC34|GC35) ModelSelected="FGC3040_34";;
	GC36|GC37) ModelSelected="FGC3040_36";;
	GC38|GC39) ModelSelected="FGC3040_38";;
	GC40|GC41) ModelSelected="FGC3545_39";;
	GC42|GC43) ModelSelected="FGC3545_41";;
	GC44|GC45) ModelSelected="FGC3545_43";;
	GC46|GC47) ModelSelected="FGC4050_46";;
	GC48|GC49) ModelSelected="FGC4050_48";;
	GC50|GC51) ModelSelected="FGC4050_50";;
	GC52|GC53) ModelSelected="FGC4560_51";;
	GC54|GC55) ModelSelected="FGC4560_53";;
	GC56|GC57) ModelSelected="FGC4560_55";;
	GC58|GC59) ModelSelected="FGC4560_57";;
	GC60|GC61) ModelSelected="FGC4560_59";;
	GC62|GC63) ModelSelected="FGC5570_61";;
	GC64|GC65) ModelSelected="FGC5570_63";;
	GC66|GC67) ModelSelected="FGC5570_65";;
	GC68|GC69) ModelSelected="FGC5570_67";;
	GC70|GC71|GC72|GC73|GC74|GC75|GC76|GC77|GC78|GC79|GC80) ModelSelected="FGC5570_67";;
esac

# Determine best model
# TODO modelName=`python3 $PROGRAMDIR/best_model.py -query $SESSIONDIR/query.fna`
modelName=/data/ppp/models/Parageobacillus_thermoglucosidasius.cnn_lstm_71/cnn_lstm.h5
# wrong model? ==> modelName=/data/ppp/models/Lactococcus_lactis.cnn_lstm_71/cnn_lstm.h5

ModelDir="Parageobacillus"
case $ModelSelected in
	General) ModelDir="Parageobacillus";;
	FGC3040_30) ModelDir="FGC3040_30.3";;
	FGC3040_32) ModelDir="FGC3040_32.3";;
	FGC3040_34) ModelDir="FGC3040_34.6";;
	FGC3040_36) ModelDir="FGC3040_36.6";;
	FGC3040_38) ModelDir="FGC3040_38.4";;
	FGC3040_40) ModelDir="FGC3040_40.3";;
	FGC3545_35) ModelDir="FGC3545_35.1";;
	FGC3545_37) ModelDir="FGC3545_37.6";;
	FGC3545_39) ModelDir="FGC3545_39.3";;
	FGC3545_41) ModelDir="FGC3545_41.4";;
	FGC3545_43) ModelDir="FGC3545_43.4";;
	FGC3545_55) ModelDir="FGC3545_45.6";;
	FGC4050_40) ModelDir="FGC4050_40.4";;
	FGC4050_42) ModelDir="FGC4050_42.1";;
	FGC4050_44) ModelDir="FGC4050_44.3";;
	FGC4050_46) ModelDir="FGC4050_46.2";;
	FGC4050_48) ModelDir="FGC4050_48.4";;
	FGC4050_50) ModelDir="FGC4050_50.2";;
	FGC4560_45) ModelDir="FGC4560_45.1";;
	FGC4560_47) ModelDir="FGC4560_47.1";;
	FGC4560_49) ModelDir="FGC4560_49.1";;
	FGC4560_51) ModelDir="FGC4560_51.1";;
	FGC4560_53) ModelDir="FGC4560_53.1";;
	FGC4560_55) ModelDir="FGC4560_55.1";;
	FGC4560_56) ModelDir="FGC4560_56.1";;
	FGC4560_57) ModelDir="FGC4560_57.1";;
	FGC4560_59) ModelDir="FGC4560_59.1";;
	FGC5570_55) ModelDir="FGC5570_55.1";;
	FGC5570_57) ModelDir="FGC5570_57.1";;
	FGC5570_59) ModelDir="FGC5570_59.1";;
	FGC5570_61) ModelDir="FGC5570_61.1";;
	FGC5570_63) ModelDir="FGC5570_63.1";;
	FGC5570_65) ModelDir="FGC5570_65.1";;
	FGC5570_67) ModelDir="FGC5570_67.1";;
	FGC5570_69) ModelDir="FGC5570_69.1";;
	FGC5570_71) ModelDir="FGC5570_71.1";;
	GC_30_35) ModelDir="GC_30_35_3";;
	GC_30_40) ModelDir="GC_30_40_4";;
	GC_35_40) ModelDir="GC_35_40_6";;
	GC_35_45) ModelDir="GC_35_45_0";;
	GC_40_45) ModelDir="GC_40_45";;
	GC_40_50) ModelDir="GC_40_50_14";;
	GC_45_50) ModelDir="GC_45_50";;
	GC_45_55) ModelDir="GC_45_55_6";;
	GC_50_55) ModelDir="GC_50_55";;
	GC_50_60) ModelDir="GC_50_60_4";;
	GC_55_60) ModelDir="GC_55_60";;
	GC_55_65) ModelDir="GC_55_65_6";;
	GC_60_65) ModelDir="GC_60_65";;
	GC_60_70) ModelDir="GC_60_70_3";;
	GC_65_70) ModelDir="GC_65_70_2";;
	Acinetobacter) ModelDir="Acinetobacter_9";;
	Bacillus) ModelDir="Bacillus_5";;
	Bacillus_amyloliquefaciens) ModelDir="Bacillus_amyloliquefaciens";;
	Bacteroides) ModelDir="Bacteroides";;
	Bradyrhizobium) ModelDir="Bradyrhizobium";;
	Campylobacter) ModelDir="Campylobacter";;
	Chlamydia) ModelDir="Chlamydia_4";;
	Chlamydophila) ModelDir="Chlamydophila_0";;
	Clostridioides) ModelDir="Clostridioides";;
	Corynebacterium) ModelDir="Corynebacterium";;
	Cupriavidus) ModelDir="Cupriavidus";;
	Cyanobacterium) ModelDir="Cyanobacterium";;
	Enterobacter) ModelDir="Enterobacter";;
	Enterococcus) ModelDir="Enterococcus_2";;
	Escherichia_coli) ModelDir="Escherichia_coli";;
	Helicobacter) ModelDir="Helicobacter";;
	Klebsiella) ModelDir="Klebsiella";;
	Lachnoclostridium) ModelDir="Lachnoclostridium";;
	LactisSpneu) ModelDir="LactisSpneu_2";;
	Lactococcus) ModelDir="Lactococcus";;
	Leptospira) ModelDir="Leptospira_9";;
	Mycobacterium) ModelDir="Mycobacterium";;
	Neisseria) ModelDir="Neisseria";;
	Parageobacillus) ModelDir="Parageobacillus";;
	Pseudomonas) ModelDir="Pseudomonas";;
	Rhodobacter) ModelDir="Rhodobacter";;
	Salmonella) ModelDir="Salmonella";;
	Shewanella) ModelDir="Shewanella";;
	Sinorhizobium) ModelDir="Sinorhizobium";;
	Staphylococcus) ModelDir="Staphylococcus_1";;
	Streptococcus) ModelDir="Streptococcus_1";;
	Streptomyces) ModelDir="Streptomyces";;
esac	


modelName=/data/ppp/models/$ModelDir/cnn_lstm.h5

# Call python script
python3 $PROGRAMDIR/ppp_Prediction_Only_Anne.py -sessiondir $SESSIONDIR -query $SESSIONDIR/query_combined_files/query -modelName $modelName -promlen 71 -pval 0.99 -out ppp



# tell webserver when run is finished
echo 'DONE' > sessionend

