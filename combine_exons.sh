#!/bin/bash

while read line
do
	echo $line
  pxcat -s aligned.header.${line}_* -o ${line}.fa
done < gene_names.txt


for filename in *.fa; do
    grep -B1 [acgt] $filename > reduced.${filename}
done


