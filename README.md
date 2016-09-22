#This readme demonstrates the steps used to train, test, and use a variation model

## Requirements:
 
- [Java 1.8+](http://www.oracle.com/technetwork/java/javase/downloads/index.html)
- wget or curl command or a web browser to download data and binaries (we use wget below)
- unzip command

## Step 1: download data and binaries
Create a working directory:
```sh
mkdir ~/variations
cd ~/variations
```
Download the jars and data (here we use wget):
```sh
wget http://gobyweb.apps.campagnelab.org/data/DLSV/model-training-1.0.2-snapshot-bin.jar
wget http://gobyweb.apps.campagnelab.org/data/DLSV/goby.jar
wget http://gobyweb.apps.campagnelab.org/data/DLSV/JGSUBWW-ZXLHITM-pickrell-NA19239_argonne-all-files.zip
wget http://gobyweb.apps.campagnelab.org/data/DLSV/VEYNFJY-HRXABYP-pickrell-NA19239_yale-all-files.zip
```
Unzip the sample alignments:
```sh
unzip JGSUBWW-ZXLHITM-pickrell-NA19239_argonne-all-files.zip -d ./somatic
unzip VEYNFJY-HRXABYP-pickrell-NA19239_yale-all-files.zip -d ./germline
```
Download and build the reference genome cache with goby:
```sh
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
java -jar goby.jar -m build-sequence-cache human_g1k_v37.fasta.gz
```
## Step 2: generate an .sbi file with goby
The .sbi file combines samples into a serialized version of each position, which will be used to generate training/testing examples:
```sh
java -jar goby.jar -m discover-sequence-variants germline/VEYNFJY-HRXABYP-pickrell-NA19239_yale.header somatic/JGSUBWW-ZXLHITM-pickrell-NA19239_argonne.header --format SEQUENCE_BASE_INFORMATION -o ./fullset --genome human_g1k_v37
```
## Step 3: generate the model
Plant mutations into the set:
```sh
java -cp model-training-1.0.2-SNAPSHOT-bin.jar org.campagnelab.dl.varanalysis.intermediaries.Mutator2 fullset.sbi mutset.sbi
```
Split the jar file into test, training, and validation:
```sh
java -cp model-training-1.0.2-SNAPSHOT-bin.jar org.campagnelab.dl.varanalysis.intermediaries.SplitFile -i mutset.sbi -f 0.8 -f 0.1 -f 0.1 -o "set_" -s train -s val -s test
```
## Step 4: train model with training set
To save time, we will only train for 10 epochs (10 full iterations over the whole training set) but in practice, you might leave out max-epochs and allow early stopping to interrupt training when validation set performance starts deteriorating.
```sh
java -cp model-training-1.0.2-SNAPSHOT-bin.jar org.campagnelab.dl.varanalysis.learning.TrainSomaticModel -t set_train.sbi -v set_val.sbi --max-epochs 10
```
## Step 5: test our model
Let's test our model on the test set. The AUC performance on the test set will be printed to standard out. The output of this test is stored in tests/[timestamp]/latest-test.tsv .
It can be used in conjuction with our [MetaR](http://metaR.campagnelab.org) software to generate ROC and Reliability curves.
```sh
mkdir tests
java -cp model-training-1.0.2-SNAPSHOT-bin.jar org.campagnelab.dl.varanalysis.learning.PredictMutations -i set_test.sbi --long-report -m ./models/[timestamp (eg 1473798870063)]
```
## Step 6:
In practice, you will want to use the model on data with no planted mutations, that way you can find actual variants. You can use Goby to do this. We'll download two pickrell 2010 datasets and look for variants.
```sh
wget http://dl.dropbox.com/u/357497/RRHCQKJ-discover-sequence-variants-demo-files.zip
unzip RRHCQKJ-discover-sequence-variants-demo-files -d ./practice
```
Create a covariates file to relate two samples to one another
```sh
echo -e "sample-id\tpatient-id\tgender\ttype\tkind-of-sample\ttissue\tparents\nZPVIZVB-pickrellNA18486_argonne\tP1\tMale\tPatient\tGermline\tBlood\tN/A\nPJCBGUJ-pickrellNA18486_yale\tP1\tMale\tPatient\tSomatic\tBlood\tN/A" > ./practice/covariates.txt
```
Now produce a vcf file with all positions obtaining a variant score (note: not true variant probability) for each position.
We will use a model threshold of .99 but you are free to use a lower threshold to investigate more possible variant sites.
```sh
java -jar goby.jar -m discover-sequence-variants practice/ZPVIZVB-pickrellNA18486_argonne.entries practice/PJCBGUJ-pickrellNA18486_yale.entries --format SOMATIC_VARIATIONS -o practice/NA18486_variants.vcf --genome human_g1k_v37 --covariates practice/covariates.txt -x SomaticVariationOutputFormat:model-path="models/[timestamp (eg 1473798870063)]/latestModel.bin" -x SomaticVariationOutputFormat:model-p-mutated-threshold:0.99
```
Our predictions will be outputted to practice/NA18486_variants.vcf.