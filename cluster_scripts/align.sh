#!/bin/bash

############      SLURM CONFIGURATION      ###################

#SBATCH --job-name=align_magnoliids_reduced
#SBATCH --partition=global
#SBATCH --account=global
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=20
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
#### 0 preparation of files and transfer to cluster
###################################################

#### Make variables with paths to directories
#Always use these paths when writing paths to files otherwise it will try to take them from your home directory

# Input directory must contain ONLY the HybPiper output fastas
# path should be "/home/ACCOUNT/HYBPIPER_OUT_FOLDER/retrieved_supercontigs/oneline/header"
# past in to path_to_dir_in
# Ensure formatted is correct (supercontigs introduce exon names into headers which can cause problems downstream)
path_to_dir_in="/data3/projects/AFRODYN2/magnoliids/hybpiper_magnoliids/retrieved_supercontigs/oneline/header_remove10";

# change output folder name
path_to_dir_out="/data3/projects/AFRODYN2/magnoliids/align_magnoliids_reduced_$SLURM_JOB_ID/";

path_to_tmp="/scratch/helmstetter_$SLURM_JOB_ID/";
	
#### Create folders on node 
echo "Transfering files to node";
mkdir $path_to_tmp
scp -r /$path_to_dir_in/* $path_to_tmp
echo "done copying files";

####
#load modules
####

module load bioinfo/mafft/7.305
module load bioinfo/Gblocks

###################################################
#### 1 Align fastas using MAFFT
###################################################

echo "starting alignment";

cd $path_to_tmp

#makes commands for all of the files in the folder and runs them in batches of 12 jobs
ls -1 ./ | \
while read sample; do
  	echo "mafft --auto ${sample} > aligned.${sample}"
done | parallel -j20

echo "done alignment";

###################################################
#### 2 Trim using Gblocks
###################################################

echo "starting trimming";

#makes commands for all of the files in the folder and runs them in batches of 12 jobs
ls -1 ./ | \
while read sample; do
	#the -b2 parameter can be removed (default) or set to 0 depending on severity of cleaning desired
  	echo "Gblocks aligned.${sample} -t=d -b5=a"
done | parallel -j20

#move to home directory
cd ~

echo "done trimming";

#Transfer output data

echo "Transfering data from node to master";

mkdir $path_to_dir_out

scp -rp $path_to_tmp/ nas3:/$path_to_dir_out/

echo "done transfer";

#### Remove tmp data on node

echo "Deleting data on node";
rm -rf $path_to_tmp
echo "Done deleting, FINISHED!";
