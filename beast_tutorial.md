# Molecular dating with BEAST

## Data preparation

### Alignments

Alignments of 75_75 loci are made as in Hybpiper pipeline but paralogs must removed and missing taxa must be filled (steps 5.1-5.3 https://github.com/ajhelmstetter/afrodyn ).

### Genetrees

Genetrees are built in Hybpiper pipeline using the genetrees.sh script. Make sure that paralog tree files are removed beforehand.

## Calculating root-to-tip variance

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

## Remove unwanted samples from alignments

If you want to remove particular individuals from a fasta alignment, put their name between the "/"

```bash
sed -i '/I02_T54/,+1 d' *.FNA
sed -i '/I02_T55/,+1 d' *.FNA
```

## Convert FASTA to nexus

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

## Model testing

Download modeltest-ng: https://github.com/ddarriba/modeltest

Open model_test.sh and edit the filepath to where you've stored modeltest-ng.
Run model_test.sh while in folder with fasta alignments

```bash
#!/bin/bash

ls -1 ./ | \
while read sample; \
do \
	#change filepath
	~/programs/modeltest-ng -d nt -i $sample -h ugif -s 3 
done
```

If you get an error saying don't have permission, run the following to make the file executable:

```bash
chmod u+x modeltest-ng
```

Then run the following to find best models as selected by BIC:

```bash
grep -A2 'Best model according to BIC' *.log | grep 'Model' 
```

Copy results into excel, add a column with numbers 0-31 and sort by name of selected model.


## Prepare XML using beauti

Download BEAST2: http://www.beast2.org/

Open beauti

file > import alignments and select your 32 chosen alignments. These should have the same number of individuals.

MAKE SURE ALIGNMENTS ARE IMPORTED IN NUMERICAL ORDER OR SUBSTITUTION MODELS WONT MATCH

Leave substitution and clock model unlinked, link trees so that we can produce a single summary tree as an output.

### Set substitution model

Use the dropdown menu to select model.

If +G set Gamma category count to 4 and make sure shape is estimated

If +I set Proportion invariant to 0.5 and check estimate

(if +I+G do both)


For models that aren't in the dropdown e.g. F81 :

https://justinbagley.rbind.io/2016/10/11/setting-dna-substitution-models-beast/


Always use empirical frequencies unless model requires otherwise (e.g. SYM)

If multiple partition ahve the same model you set the model once and then copy that model for other partitions.

After defining the model in one partition, Ctrl + click to select multiple models and clone from previously set model.

### Set clock model

Set first partition to relaxed clock log normal 

Set clock rate to 0.001

Check the box to estimate clock rate

If estimate is unchecked:

Make sure mode > automatically set clock rate is UNCHECKED

Highlight all other loci and clone from first partition

### Set priors

Tree prior = Yule model for interspecific analyses

Add prior (at the bottm) for calibration

If prior is secondary calibration on root, select all taxa then name prior root.

For secondary calibrations we might use a normal or uniform prior. For a Normal prior set mean to node age of divergence time in source tree and change Sigma value so that 2.5/97.5% covers the 95% HPD node ages.

### Set logs

Sample 10,000 trees in total.

If chain length is 10,000,000 store every 1,000. Do some test runs at this length, if this goes well then start a run  with 100,000,000 chain length sampling every 10,000.

(can find and replace in XML if many partitions)

### Run beast with multiple threads

```bash
module load bioinfo/BEAST

beast -beagle -threads 4 patch_reduced.xml
```

## STARBEAST

### Prerequisites

DAPC analysis clusters are known for each individual

* Randomly select 5 individuals per cluster using this part of the DAPC script

```R
rand5<-do.call( rbind, lapply( split(dapc_grps, dapc_grps$bar.value) ,
                               function(df) df[sample(nrow(df), min(5, nrow(df))) , ] )
)

write.csv(rand5,"dapc_grps_rand5.csv")
```

* Identify all var32 alignments and copy into a directory

* Remove individuals not part of the subset (not in the dapc_grps_rand5.csv table)

This script can be used to delete individuals from alignments:

```bash
sed -i '/I08_T55/,+1 d' *.FNA
sed -i '/I01_T10/,+1 d' *.FNA
sed -i '/I01_T12/,+1 d' *.FNA
sed -i '/I01_T13/,+1 d' *.FNA
```
Create sed commands all extra individuals and put them in a file "sed.txt

```bash
#run in directory with .FNA alignments

bash sed.txt
```

You should end up with all alignments containing (N\*5 + outgroups)\*2 number of lines.

e.g. if I have 3 clusters with 2 outgroup taxa each .FNA file should have (3\*5+2)\*2 = 34 lines

this can be checked by running 

```bash
#run in directory with .FNA alignments

wc -l *
```


* Run fasta_to_nexus.sh and model_test.sh

### load template

Open beauti and set the template to:

file > template > STARBEAST_WITH_STACEY_OPs

This improves mixing. (may need to install packages STARBEAST2)

### Import alignments

As in previous BEAST runs

###  Taxon sets

Here you assign individuals to their population genetic clusters, as inferred using DAPC,

In the species/population column enter a,b,c for clusters 1,2,3 etc.

Be sure to assign your outgroups to a different cluster "o"

### Site models 

Site models can be selected and set as in the previous BEAST tutorial

### Clock models

Clock models should be set to Strict clock for all partitions

make sure the estimate box is checked and set starting value to 0.001

### Do not touch settings in the multispecies coalesent tab

### Priors

Change the tree prior to coalescent exponential population

add a uniform prior to the root (include all taxa).

Examine your backbone tree and fine the node that represents your focal species and its sister species.

Use the 95% HPD node heights as the minimum and maximum values in the uniform prior.

### MCMC

run the chain for 300,000,000 generations, sampling tracelog, speciesTreeLogger, screenlog every 30,000 generations.

treelogs can also be sampled every 30,000 generations though it may be easier to save the xml open in sublime and find/replace "5000" with "30000"

### Launch script

You may need to install packages on your version of BEAST on the cluster for this to run. Come see me if you get errors saying "package not found" or similar.