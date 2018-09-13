Script from https://bitbucket.org/yangya/phylogenomic_dataset_construction
this script adds species/individual names to a tree file.

It can also be used on fasta files (delete hastags, but havn't tried this yet).

You need all three extra python scripts in same folder (newick3.py; phylo3.py and tree_utils.py)

The taxon name matcher file should be: 
code taxonname (tab seperated)


execute script: python2.7 taxon_name_subst.py [matchingtable] [tree file]
