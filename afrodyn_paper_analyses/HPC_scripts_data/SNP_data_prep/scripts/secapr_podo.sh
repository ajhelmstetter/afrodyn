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
#$ -N secapr_podo
############################################################
# secapr pipeline test

###################################################
#### 0 Preparation of files and transfer to cluster
###################################################

#### Creation des chemins de repertoires
# Always use these paths when writing paths to files otherwise it will try to take them from your home directory
# change depending on who is using the script/what data is being used
path_to_dir_in="/home/helmstetter/data/secapr";
path_to_scripts="/home/helmstetter/scripts";
path_to_dir_out="/data3/projects/AFRODYN2/secapr_podo_$JOB_ID/";
path_to_tmp="/scratch/helmstetter_$JOB_ID";

#### Creation du repertoire temporaire sur noeud

echo "copying files";

#make temporary directory to store files/run analyses in
mkdir $path_to_tmp

#Files that are in $path_to_dir_in:
#Annonaceae_nuc_exons.fa
#Annonaceae_pep_exons.pep
#namelist.txt #### namelist contains the sample names that will be analysed

scp nas:/$path_to_dir_in/* $path_to_tmp   ############ Modifier nas/nas2

#These steps copy all files in RUN X that have INDEX X
#Change appropriately for the files that are required
echo "copying fastqs";

scp nas:/data2/projects/raphia/RUN46/paired/trimtfiltcutRUN46_INDEX08* $path_to_tmp
scp nas:/data2/projects/raphia/RUN49/paired/trimtfiltcutRUN49_INDEX01* $path_to_tmp
scp nas:/data2/projects/raphia/catRUN46-49-50/paired/R46R49R50_TAG* $path_to_tmp
scp nas:/data2/projects/raphia/BORASSUS/RUN67/paired/trimtfiltcutRUN67-TAG-3* $path_to_tmp

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
echo "gunzipping files"
gunzip *fastq.gz
echo "done gunzipping files"

#Renames to the format e.g.I04_T44_R1
echo "renaming files"
rename -v 'trimtfiltcutRUN46_INDEX' 'I' *
rename -v 'trimtfiltcutRUN49_INDEX' 'I' *
rename -v 'trimtfiltcutRUN67' 'R67' *
rename -v 'R46R49R50' 'I12' *
rename -v -- '-TAG-' 'T' *
rename -v -- '_TAG' 'T' *
rename -v '_paired' '' *
rename -v '_R1' '_READ1' *
rename -v '_R2' '_READ2' *
echo "done renaming files"


Return to home directory

###
Format input data
###

mv $path_to_tmp/namelist_podo.txt $path_to_tmp/namelist.txt

mkdir cleaned_reads

while read name; 
do 
	echo $name
	mkdir ./cleaned_reads/"$name"_clean
	mv "$name"_READ1.fastq ./cleaned_reads/"$name"_clean
	mv "$name"_READ2.fastq ./cleaned_reads/"$name"_clean
done < $path_to_tmp/namelist.txt

cd $path_to_tmp

####
#Load modules
####

module load bioinfo/SeCaPr/1.1.4

source activate secapr_env

sleep 10s

##################################################
#### 1 Assembly
##################################################

mkdir contigs

secapr assemble_reads --input ./cleaned_reads/ --output ./contigs/ --cores 12 --disable_stats

cd contigs

ls -1 ./ | \
while read sample; \
do 
    sed -r '/^[ACGT]{,200}$/d' ${sample} | grep --no-group-separator -B1 "^[AGCT]"  > long.${sample}
done

rename 'long.' '' *

cd ../

secapr find_target_contigs --contigs ./contigs/ --reference ./palms.exons.final.fasta --output ./target_contigs --keep-duplicates --disable_stats

secapr align_sequences --sequences ./target_contigs/extracted_target_contigs_all_samples.fasta --output ./contig_alignments/ --aligner mafft --output-format fasta --ambiguous

secapr add_missing_sequences --input ./contig_alignments/ --output ./contig_alignments_no_missing/

secapr reference_assembly --reads ./cleaned_reads/ --reference_type alignment-consensus --reference ./contig_alignments --output ./remapped_reads --min_coverage 10

#secapr locus_selection --input ./remapped_reads --output ./selected_loci --n 50

secapr phase_alleles --input ./remapped_reads/ --output ./allele_sequences --min_coverage 5

#fix allele names in fasta file
sed -i -e 's/_0 |/a0 |/g' ./allele_sequences/joined_allele_fastas.fasta
sed -i -e 's/_1 |/a1 |/g' ./allele_sequences/joined_allele_fastas.fasta

secapr align_sequences --sequences ./allele_sequences/joined_allele_fastas.fasta --output ./allele_alignments/ --aligner mafft --output-format fasta --ambiguous

secapr add_missing_sequences --input ./allele_alignments/ --output ./loci_allele_alignments_complete/


##################################################
#### 5 clean up and transfer
##################################################

source deactivate

sleep 10s

module unload bioinfo/SeCaPr/1.1.4

#Transfert des donnees du noeud vers master

echo "Transfert data node -> master";

#make output folder in home directory
mkdir $path_to_dir_out

# remove cleaned reads (raw)
rm -r $path_to_tmp/cleaned_reads
scp -rp $path_to_tmp/* $path_to_dir_out 

# could move these to nas2 instead if need to store

echo "done moving";

#### Suppression du repertoire tmp noeud

echo "Deleting data on node";
#rm -rf $path_to_tmp
echo "Done deleting, FINISHED!";



