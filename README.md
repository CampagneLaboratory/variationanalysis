#This readme demonstrates the steps used to train, test, and use a variation model

## Requirements:
 
- [Java 1.8+](http://www.oracle.com/technetwork/java/javase/downloads/index.html)
- wget or curl command or a web browser to download data and binaries (we use wget below)
- unzip command

## Step 1: download data and binaries
* Create a working directory:
```sh
mkdir ~/variations
cd ~/variations
```
* Download the jars and data (here we use wget):
```sh
wget http://gobyweb.apps.campagnelab.org/data/DLSV/model-training-1.0.2-snapshot-bin.jar
wget http://gobyweb.apps.campagnelab.org/data/DLSV/goby.jar
wget http://gobyweb.apps.campagnelab.org/data/DLSV/JGSUBWW-ZXLHITM-pickrell-NA19239_argonne-all-files.zip
wget http://gobyweb.apps.campagnelab.org/data/DLSV/VEYNFJY-HRXABYP-pickrell-NA19239_yale-all-files.zip
```
_Note that you can use a more recent version of Goby 3 if available. The version bundled with this example is a preview of Goby 3.0._

* Unzip the sample alignments:
```sh
unzip JGSUBWW-ZXLHITM-pickrell-NA19239_argonne-all-files.zip -d ./somatic
unzip VEYNFJY-HRXABYP-pickrell-NA19239_yale-all-files.zip -d ./germline
```
_Note that we provide alignments in Goby format because they are much smaller than BAM, but Goby 3.0 supports BAM files as input, so you can substitute BAM alignments (.bam/.bai) in places when we show Goby alignments (.entries/.header/.index)._

* Download and build the reference genome cache with goby:
```sh
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
java -jar goby.jar -m build-sequence-cache human_g1k_v37.fasta.gz
```
## Step 2: generate an .sbi file with goby
The .sbi file combines samples into a serialized version with detailed information at each genomic position, which will be used to generate training/testing examples:
* Create the .sbi/.sbip files with the following command:
```sh
java -jar goby.jar -m discover-sequence-variants germline/VEYNFJY-HRXABYP-pickrell-NA19239_yale.header \
     somatic/JGSUBWW-ZXLHITM-pickrell-NA19239_argonne.header \
     --format SEQUENCE_BASE_INFORMATION -o ./fullset --genome human_g1k_v37
```
This command produces a .sbi and a .sbip file. The .sbi file is binary, but the .sbip is a Java properties file with information about the number of genomic sites described.
## Step 3: generate the model
* Plant mutations into the set:
```sh
java -cp model-training-1.0.2-SNAPSHOT-bin.jar org.campagnelab.dl.varanalysis.intermediaries.Mutator2 fullset.sbi \
    mutset.sbi
```
* Split the .sbi file into training (80% of training examples), validation (10%), and test (10%) datasets:
```sh
java -cp model-training-1.0.2-SNAPSHOT-bin.jar org.campagnelab.dl.varanalysis.intermediaries.SplitFile \
    -i mutset.sbi -o "set_" -s train -f 0.8   -s val -f 0.1  -s test -f 0.1
```
The command produces three datasets, called set_train.sbi, set_val.sbi and set_test.sbi.

## Step 4: train model with the training set
_To save time for this tutorial, we will only train the model for 10 epochs (10 full iterations over the whole training set).
In practice, you would leave out max-epochs and allow early stopping to interrupt training when validation set performance
starts deteriorating. This will take longer, but produce a superior model._
* Train the model for 10 epochs:
```sh
java -cp model-training-1.0.2-SNAPSHOT-bin.jar org.campagnelab.dl.varanalysis.learning.TrainSomaticModel \
    -t set_train.sbi -v set_val.sbi --max-epochs 10
```
This command generates the trained the model in ./models/[timestamp], where timestamp is a numeric value.

As training progresses, the loss is written to standard out. The performance is measured on the validation set after each epoch
is as the area under the roc curve (AUC). Training with 10 epochs should reach an AUC about 0.969.

## Step 5: test the model
Let's test the model on the test set. The AUC performance on the test set will be printed to standard out.
Additional prediction data is written to the tests directory. This file can be used to plot ROC curves and reliability diagrams.
For instance, this can be done with  [MetaR](http://metaR.campagnelab.org) and the modeling language.

```sh
mkdir tests
java -cp model-training-1.0.2-SNAPSHOT-bin.jar org.campagnelab.dl.varanalysis.learning.PredictMutations \
    -i set_test.sbi -l best \
    --long-report -m ./models/*
```
The output of this test is stored in tests/[timestamp]/best-test.tsv. The AUC of the best model on the test set
should be comparable to the AUC on the validation set.

_The model label "best" is the label of the model obtained before
performance starts to decrease on the validation set.
 If you need the very last model trained irrespective of validation performance, use -l latest_

## Step 6: use the trained model on a real dataset to find variants
In practice, you will want to use the model on data with no planted mutations, in order to identify
 actual variants. You can use Goby to do this.

* Let's download and unpack two independent pickrell 2010 datasets:

```sh
wget http://dl.dropbox.com/u/357497/RRHCQKJ-discover-sequence-variants-demo-files.zip
unzip RRHCQKJ-discover-sequence-variants-demo-files -d ./practice
```
* Create a covariates file to describe how the samples are related:
```sh
echo -e "sample-id\tpatient-id\tgender\ttype\tkind-of-sample\ttissue\tparents\nZPVIZVB-pickrellNA18486_argonne\tP1\tMale\tPatient\tGermline\tBlood\tN/A\nPJCBGUJ-pickrellNA18486_yale\tP1\tMale\tPatient\tSomatic\tBlood\tN/A" > ./practice/covariates.txt
```
This will create the following file (columns aligned for presentation, tabs must be used):
````
sample-id                       patient-id      gender  type    kind-of-sample  tissue  parents
ZPVIZVB-pickrellNA18486_argonne P1              Male    Patient Germline        Blood   N/A
PJCBGUJ-pickrellNA18486_yale    P1              Male    Patient Somatic         Blood   N/A
````
* Now produce a vcf file with all positions obtaining a model probability
 larger than 0.99.
We  use a model threshold of .99 in this tutorial, but you can lower the threshold to
investigate less confident predictions.
```sh
MODEL_TIMESTAMP=`ls -1 models`
java -jar goby.jar -m discover-sequence-variants practice/ZPVIZVB-pickrellNA18486_argonne.entries \
 practice/PJCBGUJ-pickrellNA18486_yale.entries --format SOMATIC_VARIATIONS -o practice/NA18486_variants.vcf \
 --genome human_g1k_v37 --covariates practice/covariates.txt \
 -x SomaticVariationOutputFormat:model-path=`pwd`/models/${MODEL_TIMESTAMP}/bestModel.bin \
 -x SomaticVariationOutputFormat:model-p-mutated-threshold:0.99
```
The predictions will be written to practice/NA18486_variants.vcf. The probability estimated
by the model is written to the field INFO/.