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

 1. Download and install Goby:
Add Goby to the PATH:
```sh
export PATH=${PATH}:~/IdeaProjects/git/goby3
```

2. Download and install variationanalysis
Add Variation analysis to the PATH (we assumed you downloaded version 1.2, adjust as needed):
```sh
export PATH=${PATH}:release-dlvariation_1.2/
```

### Prepare the training set

Build the Goby genome cache:
```sh
goby 10g build-sequence-cache ucsc_hg19.fasta
```


Convert the alignment file to a sequence base information (.sbi) file:
```
 goby 10g discover-sequence-variants NA12878-sorted.bam \
    --format SEQUENCE_BASE_INFORMATION     \
    -o ./NA12878 \
    --genome ucsc_hg19 \
    -n 1 -t 1   \
    --processor realign_near_indels    \
    --call-indels false -x HTSJDKReaderImpl:force-sorted=true
```
(alternatively, with GNU parallel, do 
```
export SBI_GENOME=ucsc_hg19
parallel-genotype-sbi.sh 10g NA12878-sorted.bam
```
)

Randomize sites in the sbi file:

```sh
 randomize.sh 10g -i NA12878.sbi -o NA12878-random
```

Add true labels from the ground-truth VCF:


 1. Create a genotype map (a binary file that holds the true labels).
```sh
gzip -c -d  NA12878.vcf.gz|grep -v Child >NA12878-ok.vcf
goby 4g vcf-to-genotype-map NA12878-ok.vcf -o  NA12878-true-genotypes.map
```

 2. Add true genotypes to the site information:
```sh
add-true-genotypes.sh 1g  -m true-genotypes.map -i NA12878-random.sbi  -o NA12878-random-labeled --ref-sampling-rate 0.01
```

Split sites into training, validation and test sets:

```sh
split.sh 4g -i NA12878-random-labeled -f 0.9 -s training -f 0.05 -s validation -f 0.05 -s test -o NA12878-random-labeled-
```

Train the model:
```sh
 train-genotype.sh 10g -t NA12878-random-labeled-training.sbi \
                       -v NA12878-random-labeled-validation.sbi \
                       -r 5 --mini-batch-size 512 \
                       --feature-mapper org.campagnelab.dl.genotype.mappers.GenotypeMapperV3
```

You can monitor the performance of the model being trained in the models/[timestamp] director:
```
tail -f `ls -tr1 models/*/epoch*|tail -1`
```
Will show:
```tsv
numExamples	epoch	accuracy	alleleAcc	score	condition
1187647	0	0.964349	0.992858	0.342842	not_specified
2375294	1	0.965719	0.993168	0.322654	not_specified
3562941	2	0.966069	0.993262	0.314700	not_specified
```

Training may need to proceed for tens or hundreds of epochs to converge.
On CPU, each epoch may take a few minutes. On GPU only a few seconds.
`
Calculate performance metrics on the test set:
```sh
predict-genotypes.sh 10g -m models/1482038746102 -l bestAUC -i NA12878-random-labeled-test.sbi
```

The results are in the ```predict-statistics.tsv``` file.
