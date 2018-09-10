#!/bin/bash

############      SGE CONFIGURATION      ###################
# Ecrit les erreur dans le fichier de sortie standard 
#$ -j y 

# Shell que l'on veut utiliser 
#$ -S /bin/bash 

# Email pour suivre l'execution 
#$ -M andrew.helmstetter@ird.fr

# Type de massage que l'on recoit par mail
#    -  (b) un message au demarrage
#    -  (e) a la fin
#    -  (a)  en cas d'abandon
#$ -m e 
#$ -V

# Queue que l'on veut utiliser
#$ -q bioinfo.q

# Nom du job
#$ -N align_fam
############################################################

###################################################
#### 0 preparation of files and transfer to cluster
###################################################

#### Make variables with paths to directories
#Always use these paths when writing paths to files otherwise it will try to take them from your home directory

#Input directory must contain ONLY the HybPiper output fastas.
#Ensure formatted is correct (supercontigs introduce exon names into headers which can cause problems downstream)
path_to_dir_in="/home/helmstetter/hybpiper_fam_1493584/retrieved_supercontigs/oneline/header";
path_to_dir_out="/home/helmstetter/align_fam_$JOB_ID/";
path_to_tmp="/scratch/helmstetter_$JOB_ID/";
	

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
done | parallel -j12

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
done | parallel -j12

#move to home directory
cd ~

echo "done trimming";

#Transfer output data

echo "Transfering data from node to master";

mkdir $path_to_dir_out

scp -rp $path_to_tmp/ nas:/$path_to_dir_out/

echo "done transfer";

#### Remove tmp data on node

echo "Deleting data on node";
rm -rf $path_to_tmp
echo "Done deleting, FINISHED!";
