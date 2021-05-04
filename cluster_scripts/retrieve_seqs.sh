#!/bin/bash

###################################################
#### 2 checking % sequence recovered for each exon and mapping stats
###################################################

# Runs the get_seq_lengths.py and hybpiper_stats.py parts of the hybpiper pipeline
# seem to need to be run in the directory with reads_first.py output

#cp ~/data/namelist_magnoliids_phylo.txt ./namelist.txt
#
##this script uses the reference fasta
##so if the reference is changed it must be changed in this script as well 
#bash hybpiper_stats.sh
#
####################################################
##### 3 Run intronerate and clean up folders
####################################################
#
## runs intronerate.py and cleanup.py
#bash intronerate.sh
#
#####################################################
###### 4 Paralogs
#####################################################
#
#while read i
#do
#echo $i
#python ./paralog_investigator.py $i
#done < namelist.txt
#
#####################################################
###### 5 Retrieve sequences
#####################################################

# creates folders and runs retrieve_sequences.py for output of exons (retrieved_exons)
# introns (retrieved_introns) and supercontigs (retrieved_supercontigs) 
# moves files to relevant folders

mkdir retrieved_exons
python retrieve_sequences.py Angiosperms353_targetSequences.fasta . dna
mv *.FNA retrieved_exons/

mkdir retrieved_introns
python retrieve_sequences.py Angiosperms353_targetSequences.fasta . intron
mv *.fasta retrieved_introns/

mkdir retrieved_supercontigs
mv retrieved_introns/Angiosperms353_targetSequences.fasta ./
python retrieve_sequences.py Angiosperms353_targetSequences.fasta . supercontig
mv *.fasta retrieved_supercontigs/

#This keeps all files related to paralog warnings
#Makes some weird directories and probably could be improved, but it works
find . -name '*para*' >> para.txt

mkdir retrieved_par

module load bioinfo/parallel

rm para_parallel.txt
touch para_parallel.txt

while read i;
do
	echo "mkdir -p retrieved_par/$i && cp $i retrieved_par/$i" >> para_parallel.txt
done < para.txt

#run jobs in lots of 8
parallel -j 4 < para_parallel.txt

# returns weird warnings like "cp: ommission du répertoire" for every directory, but the copy did work for the files

