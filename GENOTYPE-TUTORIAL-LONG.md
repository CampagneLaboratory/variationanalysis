### Download Data
- Start by downloading an alignment sequenced from NA12878 (from the Geuvadis project):

```
wget http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1883/E-MTAB-1883.processed.1.zip/GM_12878_No1_110407_8_sorted.bam
```

The tools demonstrated here are intended to be used with DNA-Seq alignments.
We use this sample for this tutorial because it is an RNA-Seq alignment and therefore smaller than most DNA-Seq alignments available for download.

- Download the genome sequence corresponding to this alignment:
```
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.faigz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.dict.gz
```
Note that the build is hg19. The tools demonstrated here can be used with any build corresponding to the alignment. The only requirement is that alignment and genomes must match.

- Download true genotypes from the [Platinum genome project](http://www.illumina.com/platinumgenomes/):
Click [on this link](ftp://platgene_ro@ussd-ftp.illumina.com/2016-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz) to download with a browser, then copy the downloaded file (NA12878.vcf.gz) to the project
directory.


### Download the software

 1. Download and install Goby and variationanalysis (maven is required for installation).
 We will add their tools to the path:
 (This script assumes we are at goby version 3.3.0 and variationanalysis version 1.4.0):

```
git clone https://github.com/CampagneLaboratory/goby3.git
git clone https://github.com/CampagneLaboratory/variationanalysis.git
(cd goby3; mvn install)
(cd variationanalysis; ./build-cpu.sh)
(cd goby3/snapshot-previews/; ./prepare-preview.sh 3.3.1-SNAPSHOT 1.4.1-SNAPSHOT)
cd goby3/snapshot-previews/goby-3.3.1-SNAPSHOT
export PATH=$(pwd):$PATH ; cd .. 
cd variationanalysis/bin ; export PATH=$(pwd):$PATH ; cd ../..
```

### Prepare the training set

Build the Goby genome cache:
```sh
goby 10g build-sequence-cache ucsc.hg19.fasta.gz
```

Some tools we use below execute a configure.sh file to set some of their parameters. Here are a good set of default settings.
We set all chr2 positions to be moved to seperate test dataset, and will not train our model on these examples.

```
/bin/cat <<EOM >configure.sh
export OUTPUT_PREFIX=GM_12878_gatk_realigned_sorted
export VAL_SUFFIX=validation
export MINI_BATCH_SIZE=2048
export LEARNING_RATE="3"
EVALUATION_METRIC_NAME=score
export SBI_GENOME=ucsc.hg19
export FASTA_GENOME=ucsc.hg19.fasta
export TRAINING_OPTIONS=" --genomic-context-length 29 --label-smoothing-epsilon 0.05 --memory-cache none --num-layers 4 --regularization-rate 1.6938844983121584E-11 "
export REF_SAMPLING_RATE=0.01
export REALIGN_AROUND_INDELS=false
export INCLUDE_INDELS=true
SBI_SPLIT_OVERRIDE_DESTINATION=chr2
EOM
```



Some tools are parallel, so set memory usage parameters appropriate for your machine.
Use no more than 1 thread for every 12gb of memory on your machine.
We recommend one thread per core at a maximum.
```
echo "export SBI_NUM_THREADS=5" >> configure.sh
echo "export MEM_FOR_ASYNC=16g" >> configure.sh
source configure.sh
```

#### Realign SNPs around indels
We recommend using GATK to realign around indels, for which we provide a tool that executes chromosomes in parallel.
(Make sure to have GNU parallel installed and have downloaded a GATK jar).
If you prefer, you can also use your preferred pipeline for producing a realigned bam file instead. Models trained will be very sensitive
to the specifics of the pipeline used to train them, so make sure to reproduce the pipeline for data you want to make predictions on.

```
# GATK cannot process compressed references, so we must unzip the reference and index it again:
gunzip ucsc.hg19.fasta.gz
samtools faidx ucsc.hg19.fasta

```
```
wget https://github.com/broadinstitute/picard/releases/download/2.15.0/picard.jar
java -jar picard.jar AddOrReplaceReadGroups I=GM_12878_No1_110407_8_sorted.bam O=GM_12878_No1_110407_8_sorted_addRG.bam SO=coordinate RGID=NA12878 RGLB=NA12878 RGPL=ILLUMINA RGPU=Default RGSM=NA12878 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
gunzip ucsc.hg19.fasta.gz
gunzip ucsc.hg19.fasta.fai.gz
gunzip ucsc.hg19.dict.gz
java -jar picard.jar ReorderSam I=GM_12878_No1_110407_8_sorted_addRG.bam O=GM_12878_No1_110407_8_sorted_addRG_reorder.bam REFERENCE=ucsc.hg19.fasta VALIDATION_STRINGENCY=LENIENT
parallel-gatk-realign-filtered.sh GenomeAnalysisTK.jar 20g ${NUM_THREADS} ucsc.hg19.fasta GM_12878_No1_110407_8_sorted_addRG_reorder.bam GM_12878_gatk_realigned.bam
```


Realignment can be very slow, so you can also skip this step for a faster demo:

```
mv GM_12878_No1_110407_8_sorted.bam GM_12878_gatk_realigned.bam
```

#### Convert the BAM alignment to Goby format.
(This step is optional since Goby supports BAM, Goby or CRAM aligments as input, but can speed up repeat uses of the following steps.)

Now, we will convert the realigned bam file to a Goby alignment.
```
samtools index GM_12878_gatk_realigned.bam
parallel-bam-to-goby.sh ${MEM_FOR_ASYNC} GM_12878_gatk_realigned.bam
goby 10g sort GM_12878_gatk_realigned.entries -o GM_12878_gatk_realigned_sorted

```

#### Generate the training set.
Now we will convert the alignment file to a sequence base information (.sbi) file. The sbi file serializes important
information about individual position which we train our model on. The SBI file also stores the true label for each position.
Non-variant sites will be down-sampled, because we've had better results training on a more balanced dataset.
This file will be split into training, validation, and test sets. The sets are organized by chromosome to make the test set reproducible and comparable.
```
generate-genotype-sets-0.02.sh ${MEM_FOR_ASYNC} GM_12878_gatk_realigned_sorted.entries  NA12878.vcf.gz ucsc.hg19
```


Iterate will train a model, make predictions on the hold-out test set, and report statistics about the prediction.
```
export DATE=`date +%Y-%m-%d`
echo "export DATASET=${OUTPUT_PREFIX}-${DATE}-" >> configure.sh
iterate-genotype.sh org.campagnelab.dl.genotype.mappers.GenotypeMapperV37 0
```

