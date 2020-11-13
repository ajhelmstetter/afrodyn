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
	while read name; \
	do \
		echo $name;\
		echo $sample;\
		grep -A1 ${name} $sample >> $sample.tmp;\
	done < ~/data/hybpiper/namelist_dated_beast.txt;\
done
