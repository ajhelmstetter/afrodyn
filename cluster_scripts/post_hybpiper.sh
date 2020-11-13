#!/bin/bash

# This script must be run in the "retrieved_supercontigs" directory after hybpiper finishes
# First the script makes fasta sequences on oneline
# then it removes exon names from fasta headers

ls -1 ./ | \
while read sample; \
do 
	cat $sample | awk '{if (substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' > ${sample}.oneline
done

sed -i '/^$/d' *

mkdir oneline

mv *.oneline oneline

cd oneline

rename -v '.oneline' '' *

#if line starts with > delete everything after and including -

rename -v '.fasta' '.FNA' *

for filename in *.FNA; do
	sed -r 's/^>(.{8}).*/>\1/' $filename | sed '/^>/s/-.*//' > header.${filename}
done

mkdir header

mv header.* header
