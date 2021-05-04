#!/bin/bash

############      SLURM CONFIGURATION      ###################

#SBATCH --job-name=genetrees_magnoliids_reduced
#SBATCH --partition=global
#SBATCH --account=global
#SBATCH --cpus-per-task=2
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


# phylogeny with RAxML

###################################################
#### 0 preparation of files and transfer to cluster
###################################################

#### Create paths 

#Input directory must contain the aligned and trimmed fasta files only
#CHECK FORMAT
#no exon names in supercontig fastas
#no duplicate taxa

#recommend changing output directory depending on input dataset
path_to_dir_in="/data3/projects/AFRODYN2/magnoliids/align_magnoliids_reduced_711364/helmstetter_711364/gblocks";
path_to_dir_out="/data3/projects/AFRODYN2/magnoliids/genetrees_magnoliids_reduced_$SLURM_JOB_ID";
path_to_tmp="/scratch/helmstetter_$SLURM_JOB_ID";


#### Creation du repertoire temporaire sur noeud

echo "Transferring data to node";

mkdir $path_to_tmp
scp nas3:/$path_to_dir_in/*.FNA $path_to_tmp   ############ Modifier nas/nas2
ls $path_to_tmp

echo "done copying files";

###################################################
#### 1 RAxML
###################################################

echo "starting raxml";

module load bioinfo/RAxML/8.2.9

module load bioinfo/parallel

cd $path_to_tmp

FILES=*.FNA

#make file for running analyses in parallel ### could update this to be like align (but it works now)
touch raxml_parallel.txt
for f in $FILES
do
	echo $f
	#disables the everything is undetermined warning, some trees might not have all available data.
	echo "raxmlHPC-PTHREADS -f a -x 12345 -p 12345 -T 2 -# 100 -m GTRGAMMA -O -s $f -n $f" >> raxml_parallel.txt
done

#run jobs in lots of 12 # changed from 12 to 6, as -T 2, raxml will use 2 threads, and we reserved 12 cores (slurm options)
parallel -j 20 < raxml_parallel.txt

echo "done raxml";

##################################################
####2 clean up and transfer
##################################################


#Transfer data to master
echo "Transferring data to master";

mkdir $path_to_dir_out
scp -rp $path_to_tmp/* nas3:/$path_to_dir_out/ 

echo "done moving";

#### Delete tmp data on node
#secho "Deleting data on node";
#rm -rf $path_to_tmp
#echo "Done deleting, FINISHED!";
echo "Data NOT deleted on node, FINISHED!";
echo "Must delete the folder $path_to_tmp on $SLURM_JOB_NODELIST";
