#!/bin/bash

#uses namelist to run intronerate and cleans up sample folders

while read name; 
do 
python intronerate.py --prefix $name
python cleanup.py $name
done < namelist.txt
