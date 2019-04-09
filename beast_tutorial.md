# Molecular dating with BEAST

## Data preparation

Start with alignments of 75_75 loci, paralogs removed.

### Build genetrees

As in Hybpiper using the genetrees.sh script

### Calculating root-to-tip variance

module load bioinfo/phyx

scripts downloadable at:

https://github.com/FePhyFoFum/SortaDate/blob/master/README.md

```bash
#in directory with trees

bash ../root.sh 
```

```bash
#root.sh

#!/bin/bash

#change suffix if needed
FILES=*.tre

for f in $FILES
do
	echo $f
	#replace with names of outgroups, in order of distance from focal group (-r)

	pxrr -t $f -g I12_T93,I12_T94 -r > rooted.$f
done
```

```bash
#when in folder with rooted trees
python  ~/programs/SortaDate-master/src/get_var_length.py ./ --flend .tre --outf back_roottip.txt
```

This will give output table of root-to-tip variance / total length of tree. Sort on root-to-tip variance and choose top 32 least variable loci (most clock-like)

rooted.RAxML_bipartitions.aligned.header.DN11767_10499_Q8L7R3_supercontig.FNA	0.000922119	0.322085
rooted.RAxML_bipartitions.aligned.header.DN81922_156171_O04648_supercontig.FNA	0.000374027	0.17288

You can then download your output file, open in a spreadsheet editor and sort (ascending) based on the second column (root-to-tip variance).

Copy the names of the top 32 least variable trees and put in a list:

DN11767_10499_Q8L7R3
DN81922_156171_O04648

Put the top 32 files in a folder:

```bash
cp *DN80103_66246_O22988* var32/
cp *DN79971_143901_Q9LJA3* var32/
```

### Remove unwanted samples

```bash
sed -i '/I02_T54/,+1 d' *.FNA
sed -i '/I02_T55/,+1 d' *.FNA
```

### Convert FASTA to nexus

Use script fasta_to_nexus.sh

download PGDSpider

Open program and make in PGDSpider a .spid file from fasta to nexus 

(attach spid file)

```bash
#!/bin/bash

ls -1 *.FNA | \

while read file; do
	#change path to location of PGDspider '.jar' file and '.spid' file
		java -Xmx1024m -Xms512M -jar /home/helmstet/programs/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile $file -outputfile ${file}.nexus -spid /home/helmstet/programs/PGDSpider_2.1.1.5/fasta_to_nexus.spid 
	tac ${file}.nexus | sed "1,5d" | tac > ${file}.nexus.mod
done

mkdir nexus
mv *.nexus.mod nexus

cd nexus

rename 's/.mod//' *

cd ../

rm *.nexus

```

### Model testing

Download modeltest-ng

https://github.com/ddarriba/modeltest

Run in folder with fasta alignments

If you get an error saying don't have permission, run:

```bash
chmod u+x modeltest-ng
```

to make the file executable

```bash
#!/bin/bash

ls -1 ./ | \
while read sample; \
do \
	#change filepath
	~/programs/modeltest-ng -d nt -i $sample -h ugif -s 3 
done
```
Run the following to find models selected by BIC:

```bash
grep -A2 'Best model according to BIC' *.log | grep 'Model' 
```

Copy results into excel, add a column with numbers 0-31

Sort by name of selected model


## Prepare XML using beauti

Download BEAST 2.4.4

https://github.com/CompEvol/beast2/releases?after=v2.4.6

Open beauti



file > import alignments

MAKE SURE ALIGNMENTS ARE IMPORTED IN NUMERICAL ORDER OR SUBSTITUTION MODELS WONT MATCH


### Set substitution model

For strange models e.g. F81 :

https://justinbagley.rbind.io/2016/10/11/setting-dna-substitution-models-beast/

Always use empirical frequencies unless model requires otherwise (e.g. SYM)

Set first model

Ctrl + click to select multiple models and clone from previously set model if the same substitution mode was chosen

If +G set Gamma category count to 4 and make sure shape is estimated

If +I set Proportion invariant to 0.5 and check estimate

(if +I+G do both)

### Set clock model

Set first partition to relaxed clock log normal 

Set clock rate to 0.001

If estimate is unchecked:

Make sure mode > automatically set clock rate is UNCHECKED

Highlight all other loci and clone from first partition

### Set priors

Tree = Yule model for interspecific analyses

Add prior for calibration

If prior is secondary calibration on root, select all taxa, name prior root

Set mean to node age of divergence time in other tree and change Sigma so 2.5/97.5% covers the 95% HPD node ages 

### Set logs

Sample 10,000 trees in total.

If chain length is 10,000,000 store every 1,000

At start run analysis 100,000,000 chain length

(can find and replace in XML if many partitions)

### Run beast with multiple threads

```bash
module load bioinfo/BEAST

beast -beagle -threads 6 patch_reduced.xml
```

## STARBEAST (TO BE DONE)

file > template > starbeast
