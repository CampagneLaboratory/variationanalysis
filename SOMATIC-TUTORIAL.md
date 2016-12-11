# SOMATIC-TUTORIAL.md
This tutorial demonstrates how to train, test, and use an adaptive deep learning model to predict somatic variations. 

The model is trained and tested with model-training jar, which can be compiled from source at https://github.com/CampagneLaboratory/variationanalysis (see [compilation](https://github.com/CampagneLaboratory/variationanalysis/blob/master/compiling.md) instructions).

Data used for training are processed with [Goby](http://goby.campagnelab.org) version 3.1 or later.

In this tutorial, we used RNA-Seq data described in Pickrell et al Nature. 2010 Apr 1;464(7289):768-72. doi: 10.1038/nature08872 
and publicly available from GEO accession  [GSE19480](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19480).
## Requirements:
 
- [Java 1.8+](http://www.oracle.com/technetwork/java/javase/downloads/index.html)
- wget or curl command or a web browser to download data and binaries (we use wget below)
- unzip command
- software dependencies of [DL4J](http://deeplearning4j.org/) (e.g., native libraries necessary for model training and inference, should be provided by recent Mac and Windows OS, may require specific installations on some Linux flavors).

## Step 1: download data and binaries
* Create a working directory:
```sh
mkdir ~/variations
cd ~/variations
```
* Download Goby, DL-Variation, and data (here we use wget):
```sh
wget http://gobyweb.apps.campagnelab.org/data/DLSV/release-dlvariation_1.1.1.zip [WILL WORK AFTER FORMAL RELEASE]
wget http://gobyweb.apps.campagnelab.org/data/DLSV/goby_3.1-goby.zip             [WILL WORK AFTER FORMAL RELEASE]
wget http://gobyweb.apps.campagnelab.org/data/DLSV/NA19239-all-files.zip
```
_Note that you can use a more recent version of Goby 3 if available. The version bundled with this example is a preview of Goby 3.0._

* Unzip the alignments and goby config folders:
```sh
unzip NA19239-all-files.zip
unzip release-dlvariation_1.1.1.zip
unzip goby_3.1-goby.zip
```
_Note that we provide alignments in BAM format because most people are familiar with this format. The same commands shown below work with Goby alignment. We recommend storing and processing alignments in Goby format because this is the native format for the caller, and because alignments are much smaller when stored in Goby format than in BAM.
To use Goby alignments, you can substitute  Goby alignments (.entries/.header/.index) in places of the tutorial where the command lines use BAM alignments (.bam/.bai)._

* Download and build the reference genome cache with goby:
```sh
# (can be slow, so we use a local copy instead) wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
wget http://gobyweb.apps.campagnelab.org/data/DLSV/human_g1k_v37.fasta.gz
goby-3.1/goby 4g build-sequence-cache human_g1k_v37.fasta.gz
```
## Step 2: generate an .sbi file with goby
The .sbi file combines samples into a serialized version with detailed information at each genomic position, which will be used to generate training/testing examples:
* Create the .sbi/.sbip files with the following command:
```sh
goby-3.1/goby 4g discover-sequence-variants germline/NA19239_yale.bam somatic/NA19239_argonne.bam \
     --format SEQUENCE_BASE_INFORMATION -o ./fullset --genome human_g1k_v37
```
This command should take a few minutes (there are ~15M aligned reads) and produce .sbi and .sbip file outputs. The .sbi file is binary, but the .sbip is a Java properties file with information about the number of genomic sites described.
## Step 3: generate the model
* Plant mutations into the set:
```sh
release-dlvariation_1.1.1/bin/mutate.sh 4g -i fullset.sbi -o mutset.sbi
```
* Split the .sbi file into training (80% of training examples), validation (10%), and test (10%) datasets:
```sh
release-dlvariation_1.1.1/bin/split.sh 4g \
    -i mutset.sbi -o "set_" -s train -f 0.8   -s val -f 0.1  -s test -f 0.1
```
The command produces three datasets, called set_train.sbi, set_val.sbi and set_test.sbi.

## Step 4: train model with the training set
_To save time for this tutorial, we will only train the model for 10 epochs (10 full iterations over the whole training set).
In practice, you would leave out max-epochs and allow early stopping to interrupt training when validation set performance
starts deteriorating. This will take longer, but produce a superior model._
* Train the model for 10 epochs:
```sh
release-dlvariation_1.1.1/bin/train-somatic.sh 4g \
    -t set_train.sbi -v set_val.sbi --max-epochs 10 \
    --net-architecture org.campagnelab.dl.somatic.learning.architecture.graphs.SixDenseLayersNarrower2WithFrequencyAndBase \
    --feature-mapper org.campagnelab.dl.somatic.mappers.FeatureMapperV25
```
This command generates the trained the model in ./models/[timestamp], where timestamp is a numeric value.
When the command terminates, take a look a at the content of the model directory,
as well as the model-conditions.txt file written to the current directory.

As training progresses, performance achieved by the model are written to
./models/[timestamp]/epochs-perf-log.tsv. The performance is measured on
the validation set after each epoch using the area under the roc curve (AUC).
 Another performance file records only the best AUC obtained so far:
 ./models/[timestamp]/bestAUC-perf-log.tsv so you can monitor
  this file to easily keep track of best performance. Training with 10
  epochs should reach an AUC about 0.963.

The model-conditions.txt file will contain a line like this:
```
Tag|Results|Specified_Arguments|Default_Arguments|Classname
C8CU3V|AUC 0.9630678538485847 bestModelEpoch 5 model-time 1481391664720 score 0.3521915631587974|--feature-mapper org.campagnelab.dl.somatic.mappers.FeatureMapperV25 --max-epochs 10 --net-architecture org.campagnelab.dl.somatic.learning.architecture.graphs.SixDenseLayersNarrower2WithFrequencyAndBase --training-sets [set_train.sbi] --validation-set set_val.sbi|--auc-clip-max-observations 10000 --dropout-rate null --early-stopping-measure AUC --early-stopping-num-epochs 10 --error-enrichment false --experimental-condition not_specified --feature-mapper org.campagnelab.dl.somatic.mappers.FeatureMapperV25 --gpu-device null --ignore-cache false --learning-rate 0.1 --max-epochs 10 --memory-cache false --mini-batch-size 32 --model-conditions ./model-conditions.txt --net-architecture org.campagnelab.dl.somatic.learning.architecture.graphs.SixDenseLayersNarrower2WithFrequencyAndBase --num-errors-added 16 --num-training 2147483647 --num-validation 2147483647 --parallel false --parameter-precision FP32 --previous-model-name bestAUC --previous-model-path null --random-seed 1481391664585 --regularization-rate null --training-sets [set_train.sbi] --trio false --validate-every 1 --validation-set set_val.sbi|org.campagnelab.dl.somatic.learning.TrainModelS
```
New lines will be appended each time you train a model. This is very
helpful for bookkeeping and reproducibility because it associates the
command line arguments with the performance obtained.


## Step 5: test the model
Let's test the model on the test set. The AUC performance on the test set will be printed to standard out.
Additional prediction data is written to the tests directory. This file can be used to plot ROC curves and reliability diagrams.
For instance, this can be done with  [MetaR](http://metaR.campagnelab.org) and the modeling language.

```sh
release-dlvariation_1.1.1/bin/predict.sh 4g \
    -i set_test.sbi -l bestAUC \
    -f -m ./models/*
```
The output of this test is stored in predictions/[timestamp]-bestAUC-test-set_test.tsv.
 The AUC of the best model on the test set should be comparable to the AUC on the validation set.

Predict also writes summary statistics to the file called: predict-statistics.tsv
As will model-conditions.txt, new lines are appended to this file when you run predict.sh.
This makes it easy to collate  performance of models on different datasets.

Here's the result:
```tsv
tag     prefix  auc     [auc95  auc95]  arguments
C8CU3V  bestAUC 0.961700        0.957457        0.965943        --correctness-filter null --dataset set_test.sbi --filter-auc-observations false --filter-p-max 1.0 --filter-p-min 0.0 --gpu-device null --kind test --long-report false --minibatch-size 512 --model-conditions ./model-conditions.txt --model-name bestAUC --model-path ./models/1481391664720 --num-examples 2147483647 --predict-statistics predict-statistics.tsv --records-for-auc 50000 --to-file true
```
_The model label "best" is the label of the model obtained before performance starts to decrease on the validation set.
 If you need the very last model trained irrespective of validation performance, use -l latest_

## Step 6: use the trained model on a real dataset to find variants
In practice, you will want to use the model on data with no planted mutations, in order to identify
 actual variants. You can use Goby to do this.
In this tutorial, to keep things simple, we use the files we trained with:

* We Create a [Goby covariates file](http://campagnelab.org/software/goby/tutorials/detect-somatic-variations-in-custom-patient-contigs/) to describe how the samples are related:
```sh
echo -e "sample-id\tpatient-id\tgender\ttype\tkind-of-sample\ttissue\tparents\nNA19239_yale.bam\tP1\tMale\tPatient\tGermline\tBlood\tN/A\nNA19239_argonne.bam\tP1\tMale\tPatient\tSomatic\tBlood\tN/A" > covariates.txt
```
This will create the following file (columns aligned for presentation, tabs must be used):
````
sample-id             patient-id      gender  type    kind-of-sample  tissue  parents
NA19239_yale.bam      P1              Male    Patient Germline        Blood   N/A
NA19239_argonne.bam   P1              Male    Patient Somatic         Blood   N/A
````
* Now produce a vcf file with all positions obtaining a model probability
 larger than 0.4.
We  use a model threshold of 0.4 in this tutorial so that you can see some results in the output (none are expected at the default threshold of 0.99).
This example should not produce many somatic calls at the 0.99 threshold because the samples are from the same individual and are both
technical replicates of a lymphoblastoid cell line. No somatic
variations should be found, but the model may predict a small number if
they look like real mutations. (In the tutorial data, the model may find about
25 false positive somatic variations with >99%  model probability, not a bad
false positive rate when screening about 15 million sites.)

```sh
MODEL_TIMESTAMP=`ls -1 models`
goby-3.1/goby 4g discover-sequence-variants germline/NA19239_yale.bam somatic/NA19239_argonne.bam \
 --format SOMATIC_VARIATIONS -o NA18486_variants.vcf \
 --genome human_g1k_v37 --covariates covariates.txt \
 -x SomaticVariationOutputFormat:model-path=`pwd`/models/${MODEL_TIMESTAMP}/bestAUC-ComputationGraph.bin \
 -x SomaticVariationOutputFormat:model-p-mutated-threshold:0.4
```
The predictions will be written to NA18486_variants.vcf. The probability estimated
by the model is written to the VCF field `INFO/model-probability[NA18486_argonne]`. This indicates that the probability is for a somatic mutation in sample NA18486_argonne.

## Contact
This software is developed by members of the [Campagne laboratory](http://campagnelab.org). Please address any questions or feedback to the [Goby user forum](https://groups.google.com/forum/#!forum/goby-framework), or open an issue on GitHub.
