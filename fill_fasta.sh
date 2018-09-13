#!/bin/bash
#usage bash fill_fasta.sh path/to/namelist.txt

#remove any spaces
sed -i 's/ *//g' *.FNA

#list all alignment files in folder
#remove newlines while preserving data
#outputs files with suffix .oneline
ls -1 ./*.FNA | \
while read sample; \
do 
    cat $sample | awk '{if (substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' > ${sample}.oneline
done

#get rid of introduced stuff
sed -i '/^$/d' *

#remove old files
rm *.FNA

#rename newly formatted files to same as old files
rename 's/.oneline//' *

name=$1

#files must be in current directory with suffix FNA
rename -v '.oneline' '' *

FILES=*.FNA
for f in $FILES
do
    while read line
    do
        if grep -q $line $f; then
            :
        else
            #if the name is not there, write a header
            echo -e ">${line}" >> $f
            #count number of characters on second line of alignment (bases)
            #subtract 1 because of ...newline?        
            foo=$(awk ' NR == 2' $f | wc -c)
            bases="$(($foo-1))"
            #write gap character for each base
            for ((i=1;i<=bases;i++)); do
            echo -n "-"
            done >> $f
            #newline at end of gaps
            sed -i -e '$a\' $f
        fi
    done < $name
done
