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
#$ -N hybpiper
############################################################

###################################################
#### 0 Preparation of files and transfer to cluster
###################################################

#### Creation des chemins de repertoires
# Always use these paths when writing paths to files otherwise it will try to take them from your home directory
# change depending on who is using the script/what data is being used

#Files that are in $path_to_dir_in:
#Annonaceae_nuc_exons.fa
#Annonaceae_pep_exons.pep
#namelist.txt #### namelist contains the sample names that will be analysed
path_to_dir_in="/home/helmstetter/hybpiper/data";

#this is where the hybpiper folder is located after download
#https://github.com/mossmatters/HybPiper
path_to_hybpiper="/home/helmstetter/hybpiper/program";

#The 'scripts' folder contains scripts for running intronerate
#and generating summary stats:
#hybpiper_stats.sh
#intronerate.sh
#get_seq_lengths.sh
path_to_scripts="/home/helmstetter/hybpiper/scripts";
path_to_dir_out="/home/helmstetter/output/$JOB_ID/";
path_to_tmp="/scratch/helmstetter_$JOB_ID";

#make temporary directory to store files/run analyses in
echo "copying files";
mkdir $path_to_tmp

#copy refrence files and namelist to temp directory
scp nas:/$path_to_dir_in/* $path_to_tmp

#Copy all raw fastq files needed
#Change appropriately for the files that are required
echo "copying fastqs";

scp nas2:/data/projects/afrodyn/RUN62_HISEQ/paired/INDEX12/* $path_to_tmp
scp nas2:/data/projects/afrodyn/RUN52/paired/trimtfiltcutR52-TAG-71* $path_to_tmp
scp nas2:/data/projects/afrodyn/RUN60_HISEQ/paired/trimtfiltcutRUN60_HISEQ-INDEX02-TAG-62* $path_to_tmp

echo "done copying fastqs";

#Copy scripts that are used in the pipeline to temporary folder
echo "copying scripts";

#copy all hybpiper python scripts
scp nas:/$path_to_hybpiper/*.py $path_to_tmp

#Copy shell scripts
scp nas:/$path_to_scripts/*.sh $path_to_tmp
echo "done copying scripts";

echo "done copying all files";

#Change directory to temporary directory to gunzip and rename fastqs
cd $path_to_tmp

#Gunzip files
#takes a while if you have many
echo "gunzipping files"
gunzip *fastq.gz
echo "done gunzipping files"

#Renames to the format e.g.I04_T44_R1
#made need to change patterns depending on run used
#this rename works on IRD cluster but different linux OS may require changes
echo "renaming files"
rename -v 'trimtfiltcutRUN62_HISEQ_INDEX' 'I' *
rename -v 'trimtfiltcutRUN60_HISEQ-INDEX' 'I' *
rename -v 'trimtfiltcutR52' 'R52' *
rename -v -- '-TAG-' '_T' *
rename -v '_paired' '' *
echo "done renaming files"

#load modules
module load bioinfo/SPAdes

###################################################
#### 1 reads_first
###################################################

echo "starting reads_first";

mv $path_to_tmp/namelist_i12t85.txt $path_to_tmp/namelist.txt

#A loop that reads names of samples in namelist and runs reads_first.py part of the hybpiper pipeline
while read name; 
do python $path_to_tmp/reads_first.py -b $path_to_tmp/Annonaceae_nuc_exons.fa -r $path_to_tmp/"$name"_R1.fastq $path_to_tmp/"$name"_R2.fastq --prefix $path_to_tmp/$name --bwa
done < $path_to_tmp/namelist.txt

echo "done reads_first";

###################################################
#### 2 checking % sequence recovered for each exon and mapping stats
###################################################

# Runs the get_seq_lengths.py and hybpiper_stats.py parts of the hybpiper pipeline
# seem to need to be run in the directory with reads_first.py output
# so make sure directory is $path_to_tmp before running

cd $path_to_tmp

#this script uses the reference fasta
#so if the reference is changed it must be changed in this script as well 
bash hybpiper_stats.sh

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

# creates folders and runs retrieve_sequences.py for output of exons (retrieved_exons)
# introns (retrieved_introns) and supercontigs (retrieved_supercontigs) 
# moves files to relevant folders

mkdir retrieved_exons
python retrieve_sequences.py Annonaceae_nuc_exons.fa . dna
mv *.FNA retrieved_dna/

mkdir retrieved_introns
python retrieve_sequences.py Annonaceae_nuc_exons.fa . intron
mv *.fasta retrieved_introns/

mkdir retrieved_supercontigs
python retrieve_sequences.py Annonaceae_nuc_exons.fa . supercontig
mv *.fasta retrieved_supercontigs/

#This keeps all files related to paralog warnings
#Makes some weird directories and probably could be improved, but it works
find . -name '*para*' >> para.txt

mkdir retrieved_par
while read i
do
mkdir -p retrieved_par/$i && cp $i retrieved_par/$i
done < para.txt

##################################################
#### 6 clean up and transfer
##################################################

echo "Transfert data node -> master";

#move back to home directory
cd ~

#make output folder in home directory
mkdir $path_to_dir_out

#fastq files can take up a lot of space so you may want to remove them
#rm $path_to_tmp/*.fastq

#Copies statistics and retrieved sequences to output folder in home directory
scp -p $path_to_tmp/*.txt nas:/$path_to_dir_out/
scp -rp $path_to_tmp/retrieved_dna/ nas:/$path_to_dir_out/
scp -rp $path_to_tmp/retrieved_introns/ nas:/$path_to_dir_out/
scp -rp $path_to_tmp/retrieved_supercontigs/ nas:/$path_to_dir_out/
scp -rp $path_to_tmp/retrieved_par/ nas:/$path_to_dir_out/

#copy everything, for testing
#scp -rp $path_to_tmp/* nas:/$path_to_dir_out/

echo "done moving";

#### Delete all data on node to keep space free

echo "Deleting data on node";
rm -rf $path_to_tmp
echo "Done deleting, FINISHED!";
