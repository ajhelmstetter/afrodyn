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
#$ -N RAxML_concat_fam
############################################################

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
path_to_dir_in="/home/helmstetter/annonaceae_family/concat/b2d";
path_to_dir_out="/home/helmstetter/raxml_fam_concat_$JOB_ID/";
path_to_tmp="/scratch/helmstetter_$JOB_ID";


#### Creation du repertoire temporaire sur noeud

echo "Transferring data to node";

mkdir $path_to_tmp
scp nas:/$path_to_dir_in/* $path_to_tmp   ############ Modifier nas/nas2
ls $path_to_tmp

echo "done copying files";

###################################################
#### 1 RAxML # MUST CHANGE OPTION FILEPATHS
###################################################

echo "starting raxml";

module load bioinfo/RAxML/8.2.9

cd $path_to_tmp

raxmlHPC-PTHREADS -f a -x 12345 -p 12345 -T 8 -# 100 -m GTRGAMMA -O -q ./fam_1493865.partitions -s ./fam_1493865.fa -n fam_1493865
echo "done raxml";

##################################################
####2 clean up and transfer
##################################################


#Transfer data to master
echo "Transferring data to master";

mkdir $path_to_dir_out
scp -rp $path_to_tmp/* nas:/$path_to_dir_out/ 

echo "done moving";

#### Delete tmp data on node
echo "Deleting data on node";
rm -rf $path_to_tmp
echo "Done deleting, FINISHED!";


