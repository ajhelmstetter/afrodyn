# Molecular dating with BEAST

## Data preparation

Start with alignments of 75_75 loci, paralogs removed and filled (5.1-5.3).

### Build genetrees

As in Hybpiper using the genetrees.sh script

### Calculating root-to-tip variance

Script for calculating variance is downloadable at:

https://github.com/FePhyFoFum/SortaDate/blob/master/README.md


First we need to root the trees. To run root.sh we must:

```
module load bioinfo/phyx
```

Be sure to edit the script and the outgroups to fit the dataset.

```bash
#root.sh

#!/bin/bash

#change suffix if needed
FILES=*.tre

for f in $FILES
do
	echo $f
	#replace with names of outgroups, in order of distance from focal group (-r)
	# e.g. in this example I12_T93 is further from the ingroup (closer to the root) than I12_T94
	pxrr -t $f -g I12_T93,I12_T94 -r > rooted.$f
done
```

Run the script in the directory with the tree files

```bash
bash ../root.sh 
```

Now that trees are rooted, we can calculate variance. "./" indicates that we are in the directory with the tree files. --flend is the file suffix for your trees. --outf is the outfile name

```bash
#when in folder with rooted trees
python  ~/programs/SortaDate-master/src/get_var_length.py ./ --flend .tre --outf back_roottip.txt
```

This will give output tab delimited table of root-to-tip variance / total length of tree. Sort on root-to-tip variance and choose top 32 least variable loci (most clock-like)

```
rooted.RAxML_bipartitions.aligned.header.DN11767_10499_Q8L7R3_supercontig.FNA	0.000922119	0.322085
rooted.RAxML_bipartitions.aligned.header.DN81922_156171_O04648_supercontig.FNA	0.000374027	0.17288
```

You can then download the output file, open in a spreadsheet editor and sort (ascending) based on the second column (root-to-tip variance).

Trees that are unrooted will have "NA" in this column, these can be ignored.

Copy the names of the 32 trees with the lowest root-to-tip variation (most clocklike) and paste them into a list, only reserving the exon name stem:

```
DN11767_10499_Q8L7R3
DN81922_156171_O04648
```

Copy the corresponding 32 alignments, which contain the same exon name stem, in a new folder:

```bash
mkdir var 32/

cp *DN80103_66246_O22988* var32/
cp *DN79971_143901_Q9LJA3* var32/
```

### Remove unwanted samples from alignments

If you want to remove particular individuals from a fasta alignment, put their name between the "/"

```bash
sed -i '/I02_T54/,+1 d' *.FNA
sed -i '/I02_T55/,+1 d' *.FNA
```

### Convert FASTA to nexus

Download PGDSpider: http://www.cmpg.unibe.ch/software/PGDSpider/

BEAST uses nexus format, so we need to convert our fasta alignments.

Open and make in PGDSpider a .spid file from fasta to nexus format

Edit the filepaths in fasta_to_nexus.sh for the .jar and .spid files and run

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

Download modeltest-ng: https://github.com/ddarriba/modeltest

Run in folder with fasta alignments

If you get an error saying don't have permission, run the following to make the file executable:

```bash
chmod u+x modeltest-ng
```

Edit the filepath to where you've stored modeltest-ng and run model_test.sh

```bash
#!/bin/bash

ls -1 ./ | \
while read sample; \
do \
	#change filepath
	~/programs/modeltest-ng -d nt -i $sample -h ugif -s 3 
done
```

Then run the following to find best models as selected by BIC:

```bash
grep -A2 'Best model according to BIC' *.log | grep 'Model' 
```

Copy results into excel, add a column with numbers 0-31 and sort by name of selected model.


## Prepare XML using beauti

Download BEAST2: http://www.beast2.org/

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
