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

# phylogeny with RAxML

###################################################
#### 0 preparation of files and transfer to cluster
###################################################

#### Creation des chemins de repertoires
#Always use these paths when writing paths to files otherwise it will try to take them from your home directory
path_to_dir_in="/home/helmstetter/data/anonidium_test";
path_to_scripts="/home/helmstetter/programs/HybPiper";
path_to_dir_out="/home/helmstetter/output/$JOB_ID/";
path_to_tmp="/scratch/helmstetter_$JOB_ID/";


#### Creation du repertoire temporaire sur noeud

#echo "Tranfert donnees master -> noeud";

#mkdir $path_to_tmp

#scp nas:/$path_to_dir_in/namelist.txt $path_to_tmp
#scp nas:/$path_to_dir_in/* $path_to_tmp   ############ Modifier nas/nas2

#echo "done copying files";

#echo "copying scripts";

#scp nas:/$path_to_scripts/*.py $path_to_tmp   ############ Modifier nas/nas2

#echo "done copying scripts";

#redirect to folder make pl executible

#cd $path_to_tmp/

#chmod 755 batchRAxML.pl

#echo "done moving files and excuting .pl";

####
#load modules
####

module load bioinfo/SPAdes
module load system/python/2.7.10
###################################################
#### 1 reads_first
###################################################

python /home/helmstetter/output/1006062/get_seq_lengths.py /home/helmstetter/output/1006062/Annonaceae_nuc_exons.fa /home/helmstetter/output/1006062/namelist.txt dna > /home/helmstetter/output/1006062/test_seq_lengths.txt

##################################################
####2 clean up and transfer
##################################################


#Transfert des donnees du noeud vers master

echo "Transfert data node -> master";

mkdir $path_to_dir_out
scp -rp $path_to_tmp/* nas:/$path_to_dir_out/   ############ Modifier nas/nas2

echo "done moving";

#### Suppression du repertoire tmp noeud

echo "Deleting data on node";
rm -rf $path_to_tmp
echo "Done deleting, FINISHED!";
