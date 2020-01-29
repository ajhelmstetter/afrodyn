#Demographic inference with stairwayplot

### Calculate folded SFS

for each population, should do al filtration steps e.g. remove monomorphic sites, but do not apply a MAF

Clean VCF by removing all information but alleles (1/1 0/0 ./.)

```bash
bcftools query -f '%CHROM %POS[\t%GT]\n' anni_filtered_1nm.vcf > anni_filtered_1nm.tmp

sed 's:/:        :g' anni_filtered_1nm.tmp > anni_filtered_1nm.txt

Rscript --vanilla calcSFS.r ~/programs/easySFS/anni_no_maf_no_mono/anni_filtered_1nm.txt anni1.sfs

## Above lines have been made into:
bash calcSFS.sh ~/path/to/vcf
```

### Modified filtration

We dont want to apply a MAF as that will skew the SFS

```bash
bcftools filter -i'MQ>40 && DP>25 && QD>2' anon.vcf -O vcf > anon_filtered.vcf
bcftools view -e 'AC==0 || AC==AN' -m2 -M2 -v snps -O z -o anon_filtered_st.vcf anon_filtered.vcf 
bcftools filter -i 'AVG(FMT/DP)>10' anon_filtered_st.vcf -O vcf > anon_filtered_dp.vcf
```

### Creat blueprint file
Open and edit the example blueprint file in data/stairwayplot

View the README.pdf in programs/stairway_plot_es/ for information on each line in the blueprint file 

Open VCF and copy all lengths of regions into excel and sum them to get number of sites

### Stairway plot

```bash
java -cp stairway_plot_es/ Stairbuilder anni1.blueprint

bash anon1.blueprint.sh

```