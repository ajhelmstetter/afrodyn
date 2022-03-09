#!/bin/bash

############      SLURM CONFIGURATION      ###################
#SBATCH --job-name=hybiper_magnoliids
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=12
#SBATCH --mail-user=andrew.j.helmstetter@gmail.com
#SBATCH --mail-type=ALL
############################################################

echo "JOB CONFIGURATION"
echo "Job ID: " $SLURM_JOB_ID
echo "Name of the job: " $SLURM_JOB_NAME		
echo "List of nodes allocated to the job: " $SLURM_JOB_NODELIST
echo "Number of nodes allocated to the job: " $SLURM_JOB_NUM_NODES
echo "Number of CPU tasks in this job: " $SLURM_NTASKS
echo "Directory from which sbatch was invoked: " $SLURM_SUBMIT_DIR

###################################################
#### 0 Preparation of files and transfer to cluster
###################################################

#### Creation des chemins de repertoires
# Always use these paths when writing paths to files otherwise it will try to take them from your home directory
# change depending on who is using the script/what data is being used

#Files that are in $path_to_dir_in:
#Reference file
#Angiosperms353_targetSequences.fasta
#namelist.txt #### namelist contains the sample names that will be analysed
path_to_dir_in="/data3/projects/home/helmstetter/data";

#this is where the hybpiper folder is located after download
#https://github.com/mossmatters/HybPiper
path_to_hybpiper="/data3/projects/home/helmstetter/hybpiper/program";

#The 'scripts' folder contains scripts for running intronerate
#and generating summary stats:
#hybpiper_stats.sh OR hybpiper_stats_palm.sh if working with palms 
#intronerate.sh
#get_seq_lengths.sh
path_to_scripts="/data3/projects/home/helmstetter/hybpiper/scripts";

#output folder
path_to_dir_out="/data3/projects/AFRODYN2/magnoliids/hybpiper_magnoliids_$SLURM_JOB_ID/";

#temporary folder (intermediate files)
path_to_tmp="/scratch/helmstetter_$SLURM_JOB_ID";

#make temporary directory to store files/run analyses in

echo "copying files";
mkdir $path_to_tmp

#copy refrence files and namelist to temp directory
scp nas3:$path_to_dir_in/* $path_to_tmp

#Copy scripts that are used in the pipeline to temporary folder
echo "copying scripts";

#copy all hybpiper python scripts
scp nas3:/$path_to_hybpiper/*.py $path_to_tmp

#Copy shell scripts
scp nas3:/$path_to_scripts/*.sh $path_to_tmp

echo "done copying scripts";

#Copy all raw fastq files needed
#Change appropriately for the files that are required
echo "copying fastqs";

##############################
# INSERT PATHS TO RAW FASTQs #
##############################

#load parallel
module load bioinfo/parallel

#copy all
#scp nas3:/data3/projects/AFRODYN2/magnoliids/*fastq.gz $path_to_tmp

mv $path_to_tmp/namelist_magnoliids_full.txt $path_to_tmp/namelist.txt

#move to scratch directory
cd $path_to_tmp

rm copy_parallel.txt
touch copy_parallel.txt

while read name;
do
	echo "scp nas3:/data3/projects/AFRODYN2/magnoliids/raw_data/${name}*.gz $path_to_tmp" >> copy_parallel.txt
done < namelist.txt

#run jobs in lots of4
parallel -j 4 < copy_parallel.txt

echo "done copying fastqs";

echo "done copying all files";

# gunzip and rename fastqs

echo "gunzipping files"
ls *.gz | parallel -j 12 gunzip 
echo "done gunzipping files"

#Renames to the format e.g.I04_T44_R1
#made need to change patterns depending on run used
#this rename works on IRD cluster but different linux OS may require changes
echo "renaming files"

#rename -v '_paired' '' *

echo "done renaming files"

echo "LIST OF THE FILES:"
ls

#load modules
module load bioinfo/SPAdes

###################################################
#### 1 reads_first
###################################################

echo "starting reads_first";

rm hybpiper_parallel.txt
touch hybpiper_parallel.txt

while read name;
do
	echo "python $path_to_tmp/reads_first.py -b $path_to_tmp/Angiosperms353_targetSequences.fasta -r $path_to_tmp/${name}_R1.fastq $path_to_tmp/${name}_R2.fastq --prefix $path_to_tmp/$name --cpu 1 --bwa" >> hybpiper_parallel.txt
done < namelist.txt

# /!\ Care must be taken to not overparallelise python because too many instances of python running at the same time on the same machine causes threading issues /!\
#run jobs in lots of 8
parallel -j 12 < hybpiper_parallel.txt

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

rm para_parallel.txt
touch para_parallel.txt

while read i;
do
	echo "mkdir -p retrieved_par/$i && cp $i retrieved_par/$i" >> para_parallel.txt
done < para.txt

#run jobs in lots of 8
parallel -j 12 < para_parallel.txt

# returns weird warnings like "cp: ommission du répertoire" for every directory, but the copy did work for the files

##################################################
#### 6 clean up and transfer
##################################################

echo "Transfert data node -> master";

#move back to home directory
cd ~

#make output folder in home directory
mkdir $path_to_dir_out

#fastq files can take up a lot of space so you may want to remove them
rm $path_to_tmp/*.fastq

#Copies statistics and retrieved sequences to output folder in home directory
scp -p $path_to_tmp/*.txt nas3:/$path_to_dir_out/
#scp -rp $path_to_tmp/retrieved_exons/ nas3:/$path_to_dir_out/
#scp -rp $path_to_tmp/retrieved_introns/ nas3:/$path_to_dir_out/
#scp -rp $path_to_tmp/retrieved_supercontigs/ nas3:/$path_to_dir_out/
#scp -rp $path_to_tmp/retrieved_par/ nas3:/$path_to_dir_out/

#copy everything, for testing
scp -rp $path_to_tmp/* nas:/$path_to_dir_out/

echo "done moving";

#### Delete all data on node to keep space free

#echo "Deleting data on node";
#rm -rf $path_to_tmp
#echo "Done deleting, FINISHED!";

echo "Data NOT deleted on node, FINISHED"
echo "Must delete the folder $path_to_tmp on $SLURM_JOB_NODELIST"


