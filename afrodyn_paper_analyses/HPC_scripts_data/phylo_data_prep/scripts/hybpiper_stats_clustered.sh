#!/bin/bash

python get_seq_lengths.py Annonaceae_nuc_exons_clustered.fa namelist.txt dna > test_seq_lengths.txt

python hybpiper_stats.py test_seq_lengths.txt namelist.txt > test_stats.txt

