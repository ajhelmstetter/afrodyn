# AFRODYN

Scripts and files for AFRODYN analyses

## 1 hybpiper.sh

### Rename output files and folders

Nom du job **-N** and **path_to_dir_out** should be changed for each run of the script so you know what is contained in the output folder

```bash
# Nom du job
#$ -N hybpiper_annonaceae
############################################################
```

```bash
path_to_dir_out="/home/helmstetter/output/annonaceae_$JOB_ID/";
```

### Edit filepaths of input fastqs

Replace this seciton with the full filepaths of your input fastqs:

```bash
############################
# INSERT PATHS TO RAW FASTQs 
############################

scp nas2:/data/projects/afrodyn/RUN62_HISEQ/paired/INDEX12/* $path_to_tmp
scp nas2:/data/projects/afrodyn/RUN52/paired/trimtfiltcutR52-TAG-71* $path_to_tmp
scp nas2:/data/projects/afrodyn/RUN60_HISEQ/paired/trimtfiltcutRUN60_HISEQ-INDEX02-TAG-62* $path_to_tmp
```

### Edit renaming of input fastqs

Input fastq files need to be renamed to the pattern I01_T59_R1.fastq.gz. I01 is the index or run name, T59 is the individual tag and R1 specifies read one or read two. 

Renaming should be modified to get rid of the prefixes in your fastq filenames. An example is below

```bash
echo "renaming files"

#Start with:
#trimtfiltcutRUN60_HISEQ-INDEX02-TAG-62_R2_paired.fastq.gz

rename -v 'trimtfiltcutRUN60_HISEQ-INDEX' 'I' *
rename -v -- '-TAG-' '_T' *
rename -v '_paired' '' *
echo "done renaming files"

#End with:
#I02_T62_R2.fastq.gz
```

### Make namelist, add name to hybpiper.sh and add to data folder

The namelist tells the scripts which individuals you are analysing. This file contains a list of all of the individuals you want in your analysis as follows

```bash
I02_T7
I10_T25
I10_T13
I10_T45
I12_T79
I10_T43
I10_T27
I10_T30
I10_T32
I10_T16
I10_T23
I12_T76
I10_T37
I10_T10
I12_T90
```

Place this file in the **data** folder

I suggest changing the following line of hybpiper.sh so that namelist_GROUPNAME.txt is your own unique name such as namelist_annonaceae.txt

```bash
mv $path_to_tmp/namelist_GROUPNAME.txt $path_to_tmp/namelist.txt
```

### Make and add reference file to data folder

This file contains your reference exons to be input into hybpiper.

Reference headers must have a hyphen between Species-Exon as follows:

```bash
>Monodora-DN32919_25266_Q05728
ATGAATACAATGGAGCTCTTCCATGGCATTGCTGCTGC
>Monodora-DN32989_25502_Q8GWB4
ATGTATCACTCAAGTTTTGTTAATGAAGAAGGCATTGC
>Monodora-DN33057_25972_Q93ZW3
ATGAGTGCGGAGCCCTTTTACAAGTTGGTGAAACTTGT
```

Place this file in the **data** folder

The reference must be placed in the hybpiper.sh script as well. Open up hybpiper.sh with a text editor and find & replace "Annonaceae_nuc_exons.fa" with the name of your reference file.

### Run hybpiper.sh

This script will produce stats files for the run and four folders containing exons, introns, supercontigs and one folder with the results of the paralog checks.

### post_hybpiper.sh

## 2 paralogs.sh

## 3 align.sh

### post_align.sh

## 4 genetrees.sh

### 75_75.txt

## 5 concat.sh

### fill_fasta.sh