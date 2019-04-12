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
#$ -N map_snp_anni
############################################################

###################################################
#### 0 Preparation of files and transfer to cluster
###################################################

#### Creation des chemins de repertoires
# Always use these paths when writing paths to files otherwise it will try to take them from your home directory
# change depending on who is using the script/what data is being used
path_to_ref="/home/helmstetter/annickia/ref";
path_to_names="/home/helmstetter/data/hybpiper";
path_to_dir_out="/data3/projects/AFRODYN2/map_snp_anni_$JOB_ID/";
path_to_tmp="/scratch/helmstetter_$JOB_ID";

#### Creation du repertoire temporaire sur noeud

echo "copying files";

#make temporary directory to store files/run analyses in
mkdir $path_to_tmp

#Files that are in $path_to_dir_in:
#Annonaceae_nuc_exons.fa
#Annonaceae_pep_exons.pep
#namelist.txt #### namelist contains the sample names that will be analysed

scp nas:/$path_to_ref/* $path_to_tmp   ############ Modifier nas/nas2
scp nas:/$path_to_names/* $path_to_tmp   ############ Modifier nas/nas2

#These steps copy all files in RUN X that have INDEX X
#Change appropriately for the files that are required
echo "copying fastqs";

scp nas2:/data/projects/afrodyn/RUN60_HISEQ/paired/trimtfiltcutRUN60_HISEQ-INDEX01* $path_to_tmp
scp nas2:/data/projects/afrodyn/RUN60_HISEQ/paired/trimtfiltcutRUN60_HISEQ-INDEX02* $path_to_tmp

echo "done copying fastqs";

#Copy scripts that are used in the pipeline to temporary folder
echo "copying scripts";

#copy shell scripts
#need to have:
#get_seq_lengths.sh
#hybpiper_stats.sh
#intronerate.sh

#scp nas:/$path_to_scripts/*.sh $path_to_tmp

echo "done copying scripts";

echo "done copying all files";

#Change directory to temporary directory to gunzip and rename fastqs
cd $path_to_tmp

#Renames to the format e.g.I04_T44_R1
echo "renaming files"
rename -v 'trimtfiltcutRUN60_HISEQ-INDEX' 'I' *
rename -v -- '-TAG-' '_T' *
rename -v '_paired' '' *
echo "done renaming files"

module load bioinfo/picard-tools
module load bioinfo/gatk

###################################################
#### 1 reads_first
###################################################

echo "starting reads_first";
cd $path_to_tmp

java -Xmx4g -jar $PICARD_PATH/picard.jar CreateSequenceDictionary R= joined_fasta_library.fasta O= joined_fasta_library.dict

mv $path_to_tmp/namelist_annickia.txt $path_to_tmp/namelist.txt

#A loop that reads names of samples in namelist and runs reads_first.py part of the hybpiper pipeline
while read name; do 
	echo "bwa mem $path_to_tmp/joined_fasta_library.fasta $path_to_tmp/"$name"_R1.fastq.gz $path_to_tmp/"$name"_R2.fastq.gz > "$name".sam"
done < $path_to_tmp/namelist.txt | parallel -j12

ls *.sam | parallel -j 8 'samtools sort --output-fmt=BAM {} > {.}_sorted.bam'

while read name; do 
	echo "java -Xmx4g -jar $PICARD_PATH/picard.jar MarkDuplicates I="$name"_sorted.bam O="$name"_sorted_nodup.bam M="$name"_dup.txt REMOVE_DUPLICATES=true"
done < $path_to_tmp/namelist.txt | parallel -j12

parallel -j 8 < ./parallel_readgroups.txt

parallel -j 8 < parallel_index.txt
#
parallel -j 8 < parallel_gvcf.txt
#
gatk CombineGVCFs -R joined_fasta_library.fasta -V I01_T7.g.vcf -V I01_T8.g.vcf -V I01_T9.g.vcf -V I01_T10.g.vcf -V I01_T11.g.vcf -V I01_T12.g.vcf -V I01_T13.g.vcf -V I01_T14.g.vcf -V I01_T15.g.vcf -V I01_T16.g.vcf -V I01_T17.g.vcf -V I01_T18.g.vcf -V I01_T19.g.vcf -V I01_T20.g.vcf -V I01_T21.g.vcf -V I01_T22.g.vcf -V I01_T23.g.vcf -V I01_T24.g.vcf -V I01_T25.g.vcf -V I01_T26.g.vcf -V I01_T27.g.vcf -V I01_T28.g.vcf -V I01_T29.g.vcf -V I01_T30.g.vcf -V I01_T31.g.vcf -V I01_T32.g.vcf -V I01_T33.g.vcf -V I01_T34.g.vcf -V I01_T35.g.vcf -V I01_T36.g.vcf -V I01_T37.g.vcf -V I01_T38.g.vcf -V I01_T39.g.vcf -V I01_T40.g.vcf -V I01_T41.g.vcf -V I01_T42.g.vcf -V I01_T43.g.vcf -V I01_T44.g.vcf -V I01_T45.g.vcf -V I01_T46.g.vcf -V I01_T47.g.vcf -V I01_T48.g.vcf -V I01_T49.g.vcf -V I01_T50.g.vcf -V I01_T51.g.vcf -V I01_T52.g.vcf -V I01_T53.g.vcf -V I01_T54.g.vcf -V I01_T55.g.vcf -V I01_T56.g.vcf -V I01_T57.g.vcf -V I01_T58.g.vcf -V I01_T59.g.vcf -V I01_T60.g.vcf -V I01_T61.g.vcf -V I01_T62.g.vcf -V I02_T7.g.vcf -V I02_T8.g.vcf -V I02_T9.g.vcf -V I02_T10.g.vcf -V I02_T11.g.vcf -V I02_T12.g.vcf -V I02_T13.g.vcf -V I02_T14.g.vcf -V I02_T15.g.vcf -V I02_T16.g.vcf -V I02_T17.g.vcf -V I02_T18.g.vcf -V I02_T19.g.vcf -V I02_T20.g.vcf -V I02_T21.g.vcf -V I02_T22.g.vcf -V I02_T23.g.vcf -V I02_T24.g.vcf -V I02_T25.g.vcf -V I02_T26.g.vcf -V I02_T27.g.vcf -V I02_T28.g.vcf -V I02_T29.g.vcf -V I02_T30.g.vcf -V I02_T31.g.vcf -V I02_T32.g.vcf -V I02_T33.g.vcf -V I02_T34.g.vcf -V I02_T35.g.vcf -V I02_T36.g.vcf -V I02_T37.g.vcf -V I02_T38.g.vcf -V I02_T39.g.vcf -V I02_T40.g.vcf -V I02_T41.g.vcf -V I02_T42.g.vcf -V I02_T43.g.vcf -V I02_T44.g.vcf -V I02_T45.g.vcf -V I02_T46.g.vcf -V I02_T47.g.vcf -V I02_T48.g.vcf -V I02_T49.g.vcf -V I02_T50.g.vcf -V I02_T51.g.vcf -V I02_T52.g.vcf -V I02_T53.g.vcf -V I02_T54.g.vcf -V I02_T55.g.vcf -V I02_T56.g.vcf -V I02_T57.g.vcf -V I02_T58.g.vcf -V I02_T59.g.vcf -V I02_T60.g.vcf -V I02_T61.g.vcf -V I02_T62.g.vcf -O anni.g.vcf

gatk GenotypeGVCFs -R joined_fasta_library.fasta -V anni.g.vcf -G StandardAnnotation -new-qual -O anni.vcf 

echo "done map and sort";

pwd
cd $path_to_tmp
pwd

##################################################
#### 5 clean up and transfer
##################################################

#Transfert des donnees du noeud vers master

echo "Transfert data node -> master";

#make output folder in home directory
mkdir $path_to_dir_out

#Copies statistics and retrieved sequences to output folder in home directory

scp -rp $path_to_tmp/ $path_to_dir_out

# Originally copied all output files but these were too big/many as below:
# scp -rp $path_to_tmp/I0* nas:/$path_to_dir_out/

# could move these to nas2 instead if need to store

echo "done moving";

#### Suppression du repertoire tmp noeud

echo "Deleting data on node";
rm -rf $path_to_tmp
echo "Done deleting, FINISHED!";




