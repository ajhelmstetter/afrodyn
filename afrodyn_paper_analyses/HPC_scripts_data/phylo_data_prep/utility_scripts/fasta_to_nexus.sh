#!/bin/bash

ls -1 *.FNA | \

while read file; do
	java -Xmx1024m -Xms512M -jar ~/programs/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile $file -outputfile ${file}.nexus -spid ~/programs/PGDSpider_2.1.1.5/fasta_to_nexus.spid 
	tac ${file}.nexus | sed "1,5d" | tac > ${file}.nexus.mod
done

mkdir nexus
mv *.nexus.mod nexus

cd nexus

rename 's/.mod//' *

cd ../

rm *.nexus
