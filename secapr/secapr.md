
### secapr.sh

Run without "-pe ompi X" option or it won't work.

This can produce 50gb of data, so make sure you have enough space

Namelist naming convention is different

namelist must be without "_"
e.g.

I04T44
I04T45

replace filepaths, names, email

use your own reference file



### map_snp.sh

Now we go back to hybpiper namelist (with "_" e.g. I04_T44)

find joined_fasta_library* files in secapr_output/remapped_reads/reference_seqs

copy all to your reference folder

create files

parallel_gvcf.txt
parallel_index.txt

parallel_readgroups.txt
RGID=60 should be changed for the number of the MiSeq/HiSeq run for each sample

see templates and replace names of individuals with yours

options must be changed for:
gatk CombineGVCFs 
gatk GenotypeGVCFs

### filter vcf

this will output a vcf we need to filter

bcftools version 1.9 produces an error, so use 1.8 or below (1.3 on cluster)

```bash
#typical filtering seen in papers
/home/helmstet/programs/bcftools-1.8/bcftools filter -i'MQ>40 && DP>25 && QD>2' anon.vcf -O vcf > anon_filtered.vcf

#MAF
bcftools view -q 0.01:minor anon_filtered.vcf > anon_filtered_maf.vcf

#exclude all sites at which no alternative alleles are called for any of the samples ("AC==0"), all sites at which only alternative alleles are called ("AC==AN"), and sites at which the proportion of missing data is greater than 20% ("F_MISSING > 0.2"). 
bcftools view -e 'AC==0 || AC==AN || F_MISSING > 0.2' -m2 -M2 -v snps -O z -o anon_filtered.sub2.vcf anon_filtered_maf.vcf

#Only run this  step if you want unlinked SNPs
#ensure that no SNPs are closer to each other than a minimum distance of 100 bp
vcftools --gzvcf anon_filtered.sub2.vcf --recode --thin 100000 --out anon_filtered_final.vcf
```

### DAPC clustering

Read in data with R

```R
rm(list=ls())
library(MASS)
setwd("~/Dropbox/projects/AJH_AFRODYN/annickia/")

#GOAL: create a list where you can convert indexs to names
id_ind<-read.csv("ID_index.csv")
id_vou<-read.csv("ID_voucher.csv")

only_seq<-id_vou[id_vou$ID%in%id_ind$ID,]
sort(id_ind$ID)
id_ind[order(id_ind$ID),]
foo<-cbind(only_seq[order(only_seq$ID),],id_ind[order(id_ind$ID),])
sb<-data.frame(as.character(foo$index),as.character(foo$voucher))
colnames(sb)<-c('index','voucher')
sb$index<-as.character(sb$index) 
sb$voucher<-as.character(sb$voucher)
sb<-sb[match(sort(sb$index),sb$index),]
sb<-sb[order(sb$index),]

#Read in VCF
library(pegas)
library(vcfR)
vcf <- read.vcfR("anni_filtered.sub2.vcf",checkFile = T, convertNA = T) #read in all data
head(vcf) 
vcf

### convert to genlight
aa.genlight <- vcfR2genlight(vcf, n.cores=1)
indNames(aa.genlight) <-sb$voucher
```

Follow the tutorial:

http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf

