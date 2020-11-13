#!/bin/bash

#FILES=*.FNA
#
#for f in $FILES
#do
#	while read name; 
#	do 
#		grep -A1 "$name" >> $FILES.tmp
#	done < ~/data/anonidium_test/n
#done


ls -1 ./ | \
while read sample; \
do \
	~/programs/modeltest-ng -d nt -i $sample -h ugif -s 3 
done
