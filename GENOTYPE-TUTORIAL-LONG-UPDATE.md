### Download Data
- Start by downloading an alignment sequenced from NA12878 (from the Geuvadis project):

```
wget http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1883/E-MTAB-1883.processed.1.zip/GM_12878_No1_110407_8_sorted.bam
```

The tools demonstrated here are intended to be used with DNA-Seq alignments.
We use this sample for this tutorial because it is an RNA-Seq alignment and therefore smaller than most DNA-Seq alignments available for download.

- Download the genome sequence corresponding to this alignment.
```sh
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
```
Note that the build is hg19. The tools demonstrated here can be used with any build corresponding to the alignment. The only requirement is that alignment and genomes must match.

- Download true genotypes from the [Platinum genome project](http://www.illumina.com/platinumgenomes/):
Click [on this link](ftp://platgene_ro@ussd-ftp.illumina.com/2016-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz) to download with a browser, then copy the downloaded file (NA12878.vcf.gz) to the project
directory.
Download

Assemble a single genome fasta file:
```
gzip -c -d chromFa.tar.gz |tar -xvf -
cat *.fa >ucsc_hg19.fasta
rm *.fa
```


### Download the software

 1. Download and install Goby and variationanalysis (maven is required for installation) We will add their tools to the path:
 (This script assumes we are at versions 3.2.4 and 1.2.4 respectively):
``
git clone https://github.com/CampagneLaboratory/goby3.git
git clone https://github.com/CampagneLaboratory/variationanalysis.git
(cd goby3; mvn install)
(cd variationanalysis; ./build-cpu.sh)
(goby3/formal-releases/_prepare-release.sh 3.2.4 goby3 1.2.4)
(cd release-goby_3.2.4; export PATH=$PATH:$(pwd))
(cd variationanalysis/bin; export PATH=$PATH:$(pwd))

```

### Prepare the training set

Build the Goby genome cache:
```sh
goby 10g build-sequence-cache ucsc_hg19.fasta
```

Some tools we use below execute a configure.sh file to set some of their parameters. Here are a good set of default settings.
We set all chr2 positions to be moved to seperate test dataset, and will not train our model on these examples.

```
/bin/cat <<EOM >configure.sh
export SBI_NUM_THREADS=4
export OUTPUT_PREFIX=GM_12878_realigned
export VAL_SUFFIX=validation
export MINI_BATCH_SIZE=2048
export LEARNING_RATE="3"
EVALUATION_METRIC_NAME=score
export SBI_GENOME=ucsg_hg19
export TRAINING_OPTIONS="--decision-threshold 0.4  --genomic-context-length 29 --label-smoothing-epsilon 0.17137805 --memory-cache none --num-layers 4 --regularization-rate 1.6938844983121584E-11 "
export REF_SAMPLING_RATE=0.02
export REALIGN_AROUND_INDELS=false
export INCLUDE_INDELS=true
SBI_SPLIT_OVERRIDE_DESTINATION=chr2
EOM
```



Some tools are parallelized, so set memory usage parameters appropriate for your machine.
Use no more than 1 thread for every 12gb of memory on your machine.
```
echo "export SBI_NUM_THREADS=4" >> configure.sh
echo "export MEM_FOR_ASYNC=16g" >> configure.sh
source configure.sh
```

We recommend using GATK to realign around indels, for which we provide a tool that executes chromosomes in parallel.
(Make sure to have GNU parallel installed).
If you prefer, ou can also use your preferred pipeline for producing a realigned bam file instead. Models trained will be very sensitive
to the specifics of the pipeline used to train them, so make sure to reproduce the pipeline for data you want to make predictions on.

```
parallel-gatk-realign.sh ../dependencies/gatk/GenomeAnalysisTK.jar SBI_NUM_THREADS 12g ucsc_hg19.fasta GM_12878_No1_110407_8_sorted.bam GM_12878_realigned.bam
```



Now, we will convert the realigned bam file to Goby alignements files. BAM files are not yet supported directly.
```
parallel-bam-to-goby.sh $MEM_FOR_ASYNC GM_12878_realigned.bam
```

Now we will convert the alignment file to a sequence base information (.sbi) file. The sbi file serializes important
information about individual position which we train our model on. The SBI file also stores the true label for each position.
This file will be split into training, validation, and test sets, with positions inside randomly ordered.
```
generate-genotype-sets-0.02.sh $MEM_FOR_ASYNC GM_12878_realigned  NA12878.vcf.gz ucsc_hg19
```


Iterate will train a model, make predictions on the hold-out test set, and report statistics about the prediction.
```
iterate-genotype.sh org.campagnelab.dl.genotype.mappers.GenotypeMapperV28 0
```

