This is the short version of the tutorial. It uses scripts to automate the process. If you are interested in the details of what
the scripts do, have a look at the [long version](./GENOTYPE-TUTORIAL-LONG.md).

This tutorial assumes that you are downloading data, installing software and running tools in the 
same directory. It is a good idea to create a tutorial directory and move there now. 

### Download the software

 1. Download and install Goby:
Add Goby to the PATH:
    ```sh
    wget http://gobyweb.apps.campagnelab.org/data/DLSV/goby-3.2-SNAPSHOT.zip    
    unzip goby*.zip
    export PATH=`pwd`/goby-3.2-SNAPSHOT:${PATH}
    ```

2. Download and install variationanalysis
Add the variationanalysis tools to the PATH (we assumed you downloaded version 1.2, adjust as needed):
    ```sh
    wget http://gobyweb.apps.campagnelab.org/data/DLSV/release-dlvariation_1.2-SNAPSHOT.zip 
    unzip release-dlvariation*.zip
    export PATH=`pwd`/release-dlvariation_1.2-SNAPSHOT/bin:${PATH}
    ```

### Download the Data
- Start by downloading an alignment corresponding to chr21 of sample NA12878 (from our web site):

    ```
    wget http://gobyweb.apps.campagnelab.org/data/DLSV/chr21_NA12878.zip
    unzip chr21_NA12878.zip
    ```
     _(We provide this alignment in the Goby format because it is faster to download. The software also 
    works with BAM alignments.)_
 - Download the genome sequence corresponding to this alignment.
    ```sh
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
    ```
Note that the build is hg19. The tools demonstrated here can be used with any build corresponding to the alignment.
The only requirement is that alignment and genomes must match.

 - Download true genotypes from the [Platinum genome project](http://www.illumina.com/platinumgenomes/):
Click [on this link](ftp://platgene_ro@ussd-ftp.illumina.com/2016-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz) 
to download with a browser, then copy the downloaded file (NA12878.vcf.gz) to the tutorial
directory.
on Mac:
    ````
    ftp ftp://platgene_ro@ussd-ftp.illumina.com/2016-1.0/hg19/small_variants/NA12878/NA12878.vcf.gz
    ````
    Press return when prompted for a password.
    
    Clean up the file a bit:
    ````
    gzip -c -d  NA12878.vcf.gz|grep -v Child >NA12878-ok.vcf
    ````
    
 - Assemble a single genome fasta file:
    ```
    gzip -c -d chromFa.tar.gz |tar -xvf -
    cat *.fa >ucsc_hg19.fasta
    rm *.fa
    ```

 - Index the genome
   
   Build the Goby genome cache:
   ```sh
   goby 10g build-sequence-cache ucsc_hg19.fasta
   ```
The next two steps are time-consuming, so you can skip them if you prefer and 
download a [pre-built dataset](http://gobyweb.apps.campagnelab.org/data/DLSV/chr21-NA12878-sbi-dataset.zip) 
(unzip the archive in the tutorial directory).

### Convert the VCF to a varmap
````
goby 4g vcf-to-genotype-map NA12878-ok.vcf -o  NA12878-true-genotypes.varmap
````
### Generate the training set

   ```
   export SBI_GENOME=ucsc_hg19
   export SBI_NUM_THREADS=6
   export SBI_GENOTYPE_VARMAP=NA12878-true-genotypes.varmap
   parallel-genotype-sbi.sh 10g NA12878_S1_21_dec19.entries
   ```
   _(The above assumes you can run the tool with 6 threads in parallel. Adjust if your computer
     has less cores.)_
     
### Train and evaluate the model:
In the code fragment below, adjust the DATASET variable to match the filename prefix 
for the dataset you just built (it is preset for the pre-built dataset download):
   ```sh
     export CONSIDER_INDELS=false
     export DATASET=NA12878_S1_21_dec19-2016-12-29-
     export MINI_BATCH_SIZE=2048
     iterate-genotype.sh org.campagnelab.dl.genotype.mappers.GenotypeMapperV13 1
   ```
This tool will train the model with the specified feature mapper using early stopping. 
Performance metrics will be written to the console while training proceeds. When training ends,
the model performance is evaluated on an independent test set and also reported to the console.
More data and model files are available in the model directories produced as part of training:
 
```
ls -ltr models
```

### Use the model call genotypes
The model can be used with  Goby to call genotypes and produce a VCF for new samples:

```sh
MODEL_TIMESTAMP=`ls -1 models`
goby 4g discover-sequence-variants NA12878_S1_21_dec19.entries \
 --format GENOTYPES -o NA12878_S1_chr21.vcf \
 --genome ucsc_hg19 \
 -x SomaticVariationOutputFormat:model-path=`pwd`/models/${MODEL_TIMESTAMP}/bestAUC-ComputationGraph.bin \
```
The predictions will be written to NA12878_S1_chr21.vcf. 
Note that this model is trained with too few examples to be reliable, but training with more data 
should yield reasonable models. See our preprint for details (Torracinta and Campagne, 2016. 
Training Genotype Callers with Neural Networks, BioRxiv, doi pending).