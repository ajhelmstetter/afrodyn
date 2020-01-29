# Calling SNPs and clustering with SeCaPr

## Build reference based on dataset

*secapr.sh* must be run without "-pe ompi X" option or it won't work.

This script can produce >50gb of data, so make sure you have enough space

Replace filepaths, names, email in the script and make sure the reference file is correct.

### Generate namelist

Naming convention is different to hybpiper so namelist must be without "\_"
e.g.

```
I04T44
I04T45
```

## Calling SNPs

We need to find and extract the generate reference from the *remapped_reads* folder in the secapr output and copy all associated files to our reference folder

```***
cp joined_fasta_library* ~/secapr/ref
```
To parallelize our SNP calling and VCF creation we create the following files using the templates provided:

In each file, individual names should be replaced those of your dataset

### Index SAM files
parallel_index.txt

### Modify readgroups
parallel_readgroups.txt

RGID=60 should be changed for the number of the MiSeq/HiSeq run for each sample e.g. if you have samples from run50 and run51 (R50_T1 & R51_T1) their RGID should be 50 and 51 respectively.

### Create VCF
parallel_gvcf.txt

gatk CombineGVCFs: add the names of all individual XXX.g.vcf files and the name of the combined g.vcf you want to output

gatk GenotypeGVCFs: change input g.vcf name to correspond to what you created in gatk CombineGVCFs and then a name for your final output vcf.


### Run secapr.sh

Have a look inside the script and change the relevant paths and names

Note: now we go back to *hybpiper* namelist (with "\_" e.g. I04\_T44)

### filter vcf

*secapr.sh* will output a vcf we need to filter to get high quality SNPs

bcftools version 1.9 produces an error for the first filtering step, so use 1.8 or below (1.3 on cluster)

Then reload version 1.9 for the next steps.

```bash
module load bioinfo/bcftools/1.3

#typical filtering by quality, min quality across individuals, total depth, quality by depth
bcftools filter -i 'MIN(FMT/DP)>10 && MQ>40 && DP>25 && QD>2' anni.vcf -O vcf > anni_filtered.vcf

module unload bioinfo/bcftools/1.3

module load bioinfo/bcftools/1.9 

#Minor allele filter
bcftools view -q 0.01:minor anni_filtered.vcf > anni_filtered_maf.vcf

#exclude all sites at which no alternative alleles are called for any of the samples ("AC==0"), all sites at which only alternative alleles are called ("AC==AN"), and sites at which the proportion of missing data is greater than 20% ("F_MISSING > 0.2"). 
bcftools view -e 'AC==0 || AC==AN' -m2 -M2 -v snps -O z -o anni_filtered.sub2.vcf.gz anni_filtered_maf.vcf

#Only run this step if you want unlinked SNPs (first SNP)
#ensure that no SNPs are closer to each other than a minimum distance of 100000 bp

module load bioinfo/vcftools/0.1.13

vcftools --gzvcf anon_filtered.sub2.vcf.gz --recode --thin 100000 --out anon_filtered_final.vcf
```

### Population genetic clustering

Download your filtered vcf to your own machine and open the script *clust.R*

### Read vcfR issue

There is an issue on some windows computers with read.vcfR where it thinks the file is unreadable.

If you get this error do the following:

```R
#opens package
trace(read.vcfR, edit=TRUE)

### DELETE THESE LINES
  if(file.access(file, mode = 4) != 0){
    stop(paste("File:", file, "appears to exist but is not readable!"))
  }

#save
```

Further information on DAPC can be found in the tutorial:

http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf

