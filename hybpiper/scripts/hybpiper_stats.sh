#!/bin/bash

#REPLACE $MY_REFERENCE_FILE WITH PATH TO YOUR DATASET'S OWN REFERNCE FILE

#get seq lengths of recoveries
python get_seq_lengths.py $MY_REFERENCE_FILE namelist.txt dna > test_seq_lengths.txt

#get stats
python hybpiper_stats.py test_seq_lengths.txt namelist.txt > test_stats.txt
