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
#$ -N map_snp_green
############################################################

###################################################
#### 0 Preparation of files and transfer to cluster
###################################################

#### Creation des chemins de repertoires
# Always use these paths when writing paths to files otherwise it will try to take them from your home directory
# change depending on who is using the script/what data is being used
path_to_ref="/home/helmstetter/green/ref";
path_to_names="/home/helmstetter/data/hybpiper";
path_to_dir_out="/data3/projects/AFRODYN2/map_snp_green_$JOB_ID/";
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

scp nas2:/data/projects/afrodyn/RUN60_HISEQ/paired/trimtfiltcutRUN60_HISEQ-INDEX08* $path_to_tmp
scp nas2:/data/projects/afrodyn/RUN60_HISEQ/paired/trimtfiltcutRUN60_HISEQ-INDEX09* $path_to_tmp
scp nas2:/data/projects/afrodyn/RUN53/paired/trimtfiltcutRUN53-TAG* $path_to_tmp

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

#Gunzip files

#Renames to the format e.g.I04_T44_R1
echo "renaming files"
rename -v 'trimtfiltcutRUN60_HISEQ-INDEX' 'I' *
rename -v 'trimtfiltcutRUN53' 'R53' *
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

mv $path_to_tmp/namelist_green.txt $path_to_tmp/namelist.txt

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
gatk CombineGVCFs -R joined_fasta_library.fasta  -V R53_T8.g.vcf -V R53_T9.g.vcf -V R53_T10.g.vcf -V R53_T11.g.vcf -V R53_T12.g.vcf -V R53_T13.g.vcf -V R53_T14.g.vcf -V R53_T15.g.vcf -V R53_T16.g.vcf -V R53_T17.g.vcf -V R53_T18.g.vcf -V R53_T19.g.vcf -V R53_T20.g.vcf -V R53_T21.g.vcf -V R53_T22.g.vcf -V R53_T23.g.vcf -V R53_T24.g.vcf -V R53_T25.g.vcf -V R53_T26.g.vcf -V R53_T27.g.vcf -V R53_T28.g.vcf -V R53_T29.g.vcf -V R53_T30.g.vcf -V R53_T31.g.vcf -V R53_T32.g.vcf -V R53_T33.g.vcf -V R53_T34.g.vcf -V R53_T35.g.vcf -V R53_T36.g.vcf -V R53_T37.g.vcf -V R53_T38.g.vcf -V R53_T39.g.vcf -V R53_T40.g.vcf -V R53_T41.g.vcf -V R53_T42.g.vcf -V R53_T43.g.vcf -V R53_T44.g.vcf -V R53_T45.g.vcf -V R53_T46.g.vcf -V R53_T47.g.vcf -V R53_T48.g.vcf -V R53_T49.g.vcf -V R53_T50.g.vcf -V R53_T51.g.vcf -V R53_T52.g.vcf -V R53_T53.g.vcf -V R53_T54.g.vcf -V R53_T56.g.vcf -V R53_T57.g.vcf -V R53_T58.g.vcf -V R53_T59.g.vcf -V R53_T60.g.vcf -V R53_T61.g.vcf -V R53_T62.g.vcf -V I08_T7.g.vcf -V I08_T8.g.vcf -V I08_T9.g.vcf -V I08_T10.g.vcf -V I08_T11.g.vcf -V I08_T12.g.vcf -V I08_T13.g.vcf -V I08_T14.g.vcf -V I08_T15.g.vcf -V I08_T16.g.vcf -V I08_T17.g.vcf -V I08_T18.g.vcf -V I08_T19.g.vcf -V I08_T20.g.vcf -V I08_T21.g.vcf -V I08_T22.g.vcf -V I08_T23.g.vcf -V I08_T24.g.vcf -V I08_T25.g.vcf -V I08_T26.g.vcf -V I08_T27.g.vcf -V I08_T28.g.vcf -V I08_T29.g.vcf -V I08_T30.g.vcf -V I08_T31.g.vcf -V I08_T32.g.vcf -V I08_T33.g.vcf -V I08_T34.g.vcf -V I08_T35.g.vcf -V I08_T36.g.vcf -V I08_T37.g.vcf -V I08_T38.g.vcf -V I08_T39.g.vcf -V I08_T40.g.vcf -V I08_T41.g.vcf -V I08_T42.g.vcf -V I08_T43.g.vcf -V I09_T7.g.vcf -V I09_T8.g.vcf -V I09_T9.g.vcf -V I09_T10.g.vcf -V I09_T11.g.vcf -V I09_T12.g.vcf -V I09_T13.g.vcf -V I09_T14.g.vcf -V I09_T15.g.vcf -V I09_T16.g.vcf -V I09_T17.g.vcf -V I09_T19.g.vcf -V I09_T20.g.vcf -V I09_T21.g.vcf -V I09_T22.g.vcf -V I09_T23.g.vcf -V I09_T24.g.vcf -V I09_T25.g.vcf -V I09_T27.g.vcf -V I09_T28.g.vcf -V I09_T29.g.vcf -V I09_T30.g.vcf -V I09_T31.g.vcf -V I09_T32.g.vcf -V I09_T33.g.vcf -V I09_T34.g.vcf -V I09_T35.g.vcf -V I09_T36.g.vcf -V I09_T37.g.vcf -V I09_T38.g.vcf -V I09_T39.g.vcf -V I09_T40.g.vcf -V I09_T41.g.vcf -V I09_T42.g.vcf -V I09_T43.g.vcf -V I09_T44.g.vcf -V I09_T45.g.vcf -V I09_T46.g.vcf -V I09_T47.g.vcf -V I09_T48.g.vcf -V I09_T49.g.vcf -V I09_T50.g.vcf -V I09_T51.g.vcf -V I09_T52.g.vcf -V I09_T53.g.vcf -V I09_T54.g.vcf -V I09_T55.g.vcf -V I09_T56.g.vcf -V I09_T57.g.vcf -V I09_T58.g.vcf -V I09_T59.g.vcf -V I09_T60.g.vcf -V I09_T61.g.vcf -V I09_T62.g.vcf -O green.g.vcf

gatk GenotypeGVCFs -R joined_fasta_library.fasta -V green.g.vcf -G StandardAnnotation -new-qual -O green.vcf 

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



