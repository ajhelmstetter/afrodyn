#!/bin/bash

############      SLURM CONFIGURATION      ###################

#SBATCH --job-name=paralogs_TEST
#SBATCH --partition=highmem
#SBATCH --nodelist=node4
# --account=project_group
#SBATCH --cpus-per-task=2 
#SBATCH --ntasks-per-node=4
#SBATCH --mail-user=leo-paul.dagallier@ird.fr
#SBATCH --mail-type=ALL

############################################################

# Return the JOB infos in the log file:
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
path_to_dir_in="/data3/projects/AFRODYN2/hybpiper_back_1590007/retrieved_par";
path_to_dir_out="/home/helmstetter/paralog_back_$JOB_ID/";
path_to_tmp="/scratch/helmstetter_$JOB_ID";

#### Creation du repertoire temporaire sur noeud

echo "Transferring data to node";

mkdir $path_to_tmp
scp -rp $path_to_dir_in/* $path_to_tmp   ############ Modifier nas/nas2
ls $path_to_tmp

echo "done copying files";

module load bioinfo/RAxML/8.2.9
module load bioinfo/mafft/7.305

cd $path_to_tmp

#Previously par_table.sh 
#Finds all genes with warnings from hybpiper run and makes a list 

find . -name genes_with_paralog_warnings.txt >> par_list.txt

rm para_table.txt
touch para_table.txt

while read i
do
cat $i >> para_table.txt
echo $i >> para_table.txt

echo -e "\n" >> para_table.txt
done < par_list.txt

sed -i 's/[.].*//' para_table.txt # removes every character after the "." character 
sed -i '/^$/d' para_table.txt # removes empty lines

cat para_table.txt | sort -f | uniq > loci_list.txt

# Combine fasta for each loci
echo "combine fasta for each loci";
while read i
do
	echo $i
	rm -r ${i}_fastas
	rm ${i}_paralog_fasta_list.txt
	
	find . -name ${i}_paralogs.fasta > ${i}_paralog_fasta_list.txt
	
	mkdir ${i}_fastas
	
	while read j
	do
		nom=$(echo $j | cut -c3-9)
		cp $j ${i}_fastas/${nom}.fasta
	done < ${i}_paralog_fasta_list.txt

	cat ${i}_fastas/* >> ${i}_fastas/combined_${i}.fasta
done < loci_list.txt
echo "done combining fasta for each loci";

# Gather the combined fastas into a single directory and run alignments
echo "gather combined fasta";
	mkdir combined_fastas
	find . -name 'combined_*.fasta' -exec cp -t combined_fastas {} +
echo "done gathering combined fasta";

echo "starting alignment";
	cd combined_fastas
	#makes commands for all of the files in the folder and runs them in batches of jobs
	ls -1 ./ | \
		while read sample; do
		  	echo "mafft --auto ${sample} > aligned.${sample}"
		done | parallel -j4 #change depending on dataset size: has to be the same as --ntasks-per-node (parameter in SLURM configuration)
echo "done alignment";

# Gather the aligned combined fastas into a single directory and run RAxML
echo "gather aligned combined fasta";
mkdir aligned_combined_fastas
find . -name 'aligned.combined_*' -exec cp -t aligned_combined_fastas {} +
echo "done gathering aligned combined fasta";

echo "starting raxml";
cd aligned_combined_fastas
	#makes commands for all of the files in the folder and runs them in batches of jobs
	ls -1 ./ | \
		while read sample; do
		  	echo "raxmlHPC-PTHREADS -f a -x 12345 -p 12345 -T 2 -# 100 -m GTRGAMMA -O -s ${sample} -n ${sample}" #change -T depending on dataset size: has to be the same as the --cpus-per-task parameter (in SLURM configuration)
		done | parallel -j4 #change depending on dataset size: has to be the same as the --ntasks-per-node parameter (in SLURM configuration)
echo "done raxml";

# Gather the RAxML trees into a single directory (to be transfered locally for plot_paralogs.R)
echo "gather the trees"; 
cd $path_to_tmp
mkdir trees
find . -name '*bipartitions.*' -exec cp -t trees {} +
echo "trees gathered";


#Transfer output data

echo "Transfering data from node to master";

mkdir $path_to_dir_out

scp -rp $path_to_tmp/ nas:/$path_to_dir_out/

echo "done transfer";

#### Remove tmp data on node

echo "Deleting data on node";
rm -rf $path_to_tmp
echo "Done deleting, FINISHED!";
