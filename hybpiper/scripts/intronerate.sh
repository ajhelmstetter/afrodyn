#!/bin/bash

# A for loop that iterates over samples 7 to 62 in index 3 and 4
# and runs intronerate.py and cleanup.py
# Change sequence for i and indices depending which samples need to be analysed
#
#for i in {7..70}
#do
#python intronerate.py --prefix I06_T$i
#python cleanup.py I06_T$i
#done
#
#for i in {7..62}
#do
#python intronerate.py --prefix I10_T$i
#python cleanup.py I10_T$i
#done
#
#for i in {71..102}
#do
#python intronerate.py --prefix I12_T$i
#python cleanup.py I12_T$i
#done

#uses namelist to run intronerate

while read name; 
do 
python intronerate.py --prefix $name
python cleanup.py $name
done < namelist.txt