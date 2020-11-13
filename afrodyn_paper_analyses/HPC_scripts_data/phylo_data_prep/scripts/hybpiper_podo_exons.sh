#!/bin/bash

############      SGE CONFIGURATION      ###################
# Ecrit les erreur dans le fichier de sortie standard 
#$ -j y 

# Shell que l'on veut utiliser 
#$ -S /bin/bash 

# Email pour suivre l'execution 
#$ -M andrew.j.helmstetter@gmail.com

# Type de massage que l'on recoit par mail
#    -  (b) un message au demarrage
#    -  (e) a la fin
#    -  (a)  en cas d'abandon
#$ -m e 
#$ -V

# Queue que l'on veut utiliser
#$ -q bioinfo.q

# Nom du job
#$ -N hybpiper_podo_a_exons
############################################################
# hybpiper pipeline test

###################################################
#### 0 Preparation of files and transfer to cluster
###################################################

#### Creation des chemins de repertoires
# Always use these paths when writing paths to files otherwise it will try to take them from your home directory
# change depending on who is using the script/what data is being used
path_to_dir_in="/home/helmstetter/data/hybpiper";
path_to_hybpiper="/home/helmstetter/programs/HybPiper";
path_to_scripts="/home/helmstetter/scripts";
path_to_dir_out="/home/helmstetter/hybpiper_podo_a_exons_$JOB_ID/";
path_to_tmp="/scratch/helmstetter_$JOB_ID";

#### Creation du repertoire temporaire sur noeud

echo "copying files";

#make temporary directory to store files/run analyses in
mkdir $path_to_tmp

#Files that are in $path_to_dir_in:
#Annonaceae_nuc_exons.fa
#Annonaceae_pep_exons.pep
#OR palms.exons_split.final.fastasta
#namelist.txt #### namelist contains the sample names that will be analysed

scp nas:/$path_to_dir_in/* $path_to_tmp   ############ Modifier nas/nas2

#These steps copy all files in RUN X that have INDEX X
#Change appropriately for the files that are required
echo "copying fastqs";

scp nas:/data2/projects/raphia/RUN46/paired/trimtfiltcutRUN46_INDEX08* $path_to_tmp
scp nas:/data2/projects/raphia/RUN49/paired/trimtfiltcutRUN49_INDEX01* $path_to_tmp
scp nas:/data2/projects/raphia/catRUN46-49-50/paired/R46R49R50_TAG* $path_to_tmp
scp nas:/data2/projects/raphia/BORASSUS/RUN67/paired/trimtfiltcutRUN67-TAG-3* $path_to_tmp

#acaulis
scp nas:/data2/projects/raphia/RUN49/paired/trimtfiltcutRUN49_INDEX01-TAG-22_* $path_to_tmp
scp nas:/data2/projects/raphia/RUN49/paired/trimtfiltcutRUN49_INDEX01-TAG-39_* $path_to_tmp

echo "done copying fastqs";

#Copy scripts that are used in the pipeline to temporary folder
echo "copying scripts";

#copy all hybpiper python scripts
scp nas:/$path_to_hybpiper/*.py $path_to_tmp   ############ Modifier nas/nas2

#copy shell scripts
#need to have:
#get_seq_lengths.sh
#hybpiper_stats.sh
#intronerate.sh

scp nas:/$path_to_scripts/*.sh $path_to_tmp

echo "done copying scripts";

echo "done copying all files";

#Change directory to temporary directory to gunzip and rename fastqs
cd $path_to_tmp

#Gunzip files
echo "gunzipping files"
gunzip *fastq.gz
echo "done gunzipping files"

#Renames to the format e.g.I04_T44_R1
echo "renaming files"
rename -v 'trimtfiltcutRUN46_INDEX' 'I' *
rename -v 'trimtfiltcutRUN49_INDEX' 'I' *
rename -v 'trimtfiltcutRUN67' 'R67' *
rename -v 'R46R49R50' 'I12' *
rename -v -- '-TAG-' '_T' *
rename -v '_TAG' '_T' *
rename -v '_paired' '' *
echo "done renaming files"

#Return to home directory
cd ~

####
#load modules
####

module load bioinfo/SPAdes

###################################################
#### 1 reads_first
###################################################

echo "starting reads_first";

mv $path_to_tmp/namelist_podo_a.txt $path_to_tmp/namelist.txt
mv $path_to_tmp/palms.exons_split.final.fasta $path_to_tmp/palms.exons_split.final.fasta 

#A loop that reads names of samples in namelist and runs reads_first.py part of the hybpiper pipeline
while read name; 
do python $path_to_tmp/reads_first.py -b $path_to_tmp/palms.exons_split.final.fasta -r $path_to_tmp/"$name"_R1.fastq $path_to_tmp/"$name"_R2.fastq --prefix $path_to_tmp/$name --bwa
done < $path_to_tmp/namelist.txt

echo "done reads_first";

echo "done reads_first";

###################################################
#### 2 checking % sequence recovered for each exon and mapping stats
###################################################

# Runs the get_seq_lengths.py and hybpiper_stats.py parts of the hybpiper pipeline
# seem to need to be run in the directory with reads_first.py output
# so change directory to $path_to_tmp before running

pwd
cd $path_to_tmp
pwd

bash hybpiper_stats_palm_exons.sh

###################################################
#### 3 Run intronerate and clean up folders
###################################################

# runs intronerate.py and cleanup.py

bash intronerate.sh
####################################################
##### 4 Paralogs
####################################################

while read i
do
echo $i
python ./paralog_investigator.py $i
done < namelist.txt

####################################################
##### 5 Retrieve sequences
####################################################

# creates folders and runs retrieve_sequences.py for output of exons (retrieved_dna)
# introns (retrieved_introns) and supercontigs (retrieved_supercontigs) 
# moves files to relevant folders

mkdir retrieved_dna
python retrieve_sequences.py palms.exons_split.final.fasta . dna
mv *.FNA retrieved_dna/

mkdir retrieved_introns
python retrieve_sequences.py palms.exons_split.final.fasta . intron
mv *.fasta retrieved_introns/

mkdir retrieved_supercontigs
python retrieve_sequences.py palms.exons_split.final.fasta . supercontig
mv *.fasta retrieved_supercontigs/


find . -name '*para*' >> para.txt

mkdir retrieved_par
while read i
do
mkdir -p retrieved_par/$i && cp $i retrieved_par/$i
done < para.txt

#move back to home directory
cd ~

##################################################
#### 6 clean up and transfer
##################################################

#Transfert des donnees du noeud vers master

echo "Transfert data node -> master";

#make output folder in home directory
mkdir $path_to_dir_out

#remove fastq files
#rm $path_to_tmp/*.fastq

#Copies statistics and retrieved sequences to output folder in home directory
scp -p $path_to_tmp/*.txt nas:/$path_to_dir_out/
scp -rp $path_to_tmp/retrieved_dna/ nas:/$path_to_dir_out/
scp -rp $path_to_tmp/retrieved_introns/ nas:/$path_to_dir_out/
scp -rp $path_to_tmp/retrieved_supercontigs/ nas:/$path_to_dir_out/
scp -rp $path_to_tmp/retrieved_par/ nas:/$path_to_dir_out/

#copy everything, for testing
#scp -rp $path_to_tmp/* nas:/$path_to_dir_out/

# Originally copied all output files but these were too big/many as below:
# scp -rp $path_to_tmp/I0* nas:/$path_to_dir_out/
# could move these to nas2 instead if need to store

echo "done moving";

#### Delete all data on node to keep space free

echo "Deleting data on node";
rm -rf $path_to_tmp
echo "Done deleting, FINISHED!";