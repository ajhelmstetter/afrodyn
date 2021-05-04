#!/bin/bash
#put paralogs/filtering scripts with commands in genetrees output folder
#run in genetrees output folder
#usage bash post_genetrees_astral.sh /path/to/astral.jar

#remove old trees directory
rm -r trees
mkdir trees

#find tree files (with support)
find . -name '*bipartitions.*' -exec cp -t trees {} +

cd trees

#Uncomment to run paralog removal
mkdir paralogs
bash ../mv_paralogs_handpicked.txt

#Uncomment to run ASTRAL on filtered loci
mkdir filtered
bash ../cp_filtered.txt
cd filtered/

#combine all trees
cat *.FNA > all.trees

#Collapse branches with bootstrap < 10
module load bioinfo/newick-utils/1.6

nw_ed all.trees 'i & b<=10' o > all_bs10.trees

#run ASTRAL with LPP
java -jar $1 -t 3 -i all_bs10.trees -o astral_all_bs10_LPP.tre 2> astral_all_bs10_LPP.log

#run ASTRAL with quartet scores
java -jar $1 -t 1 -i all_bs10.trees -o astral_all_bs10_QS.tre 2> astral_all_bs10_QS.log

#bring together output
mkdir $2

mv *.trees $2
mv *.tre $2
mv *.log $2

