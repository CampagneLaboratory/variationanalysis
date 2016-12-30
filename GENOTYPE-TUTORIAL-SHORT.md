This is the short version of the tutorial. It uses scripts to automate the process. If you are interested in the details of what
the scripts do, have a look at the [long version](./GENOTYPE-TUTORIAL-LONG.md).

This tutorial assumes that you are downloading data, installing software and running tools in the 
same directory. It is a good idea to create a tutorial directory and move there now. 

### Download the software

 1. Download and install Goby:
Add Goby to the PATH:
    ```sh
    wget http://gobyweb.apps.campagnelab.org/data/DLSV/goby_3.2-SNAPSHOT-bin.zip    
    unzip goby*.zip
    export PATH=${PATH}:~/IdeaProjects/git/goby3
    ```

2. Download and install variationanalysis
Add the variationanalysis tools to the PATH (we assumed you downloaded version 1.2, adjust as needed):
    ```sh
    wget wget http://gobyweb.apps.campagnelab.org/data/DLSV/release-dlvariation_1.2-SNAPSHOT.zip 
    unzip release-dlvariation*.zip
    export PATH=${PATH}:release-dlvariation-1.2-SNAPSHOT/bin
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


### Generate the training set
   ```
   export SBI_GENOME=ucsc_hg19
   export SBI_NUM_THREADS=6
   parallel-genotype-sbi.sh 10g NA12878_S1.bam
   ```
   _(The above assumes you can run the tool with 6 threads in parallel. Adjust if your computer
     has less cores.)_
     
### Train and evaluate the model:

   ```sh
     export CONSIDER_INDELS=false
     iterate-genotype.sh org.campagnelab.dl.genotype.mappers.GenotypeMapperV13 1
   ```
This tool will train the model with the specified feature mapper using early stopping. 
Performance metrics will be written to the console while training proceeds. When training ends,
the model performance is evaluated on an independent test set and also reported to the console.
More data and model files are available in the model directories produced as part of training:
 
```
ls -ltr models
```
