##CHANGE LOG

### 1.4.0 (Nov 2017)
* Various bug fixes. Indel performance is now state of the art on NA12878 with the V37 mapper 
  (org.campagnelab.dl.genotype.mappers.GenotypeMapperV37). Note that .sbi files must be 
  rebuilt with Goby 3.3.0+ (some bug fixes were done in Goby).
* Framework TrainModel: added support for Mixup (See https://arxiv.org/abs/1710.09412.)
  Training with mixup is supported for all types of models developed in this project.
* Debugged LSTM inputs for genotypes. Similar performance to fully connected architecture,
  but more compact models. Use V38 genotype mapper with the LSTM architecture. 
* Draft of genotype segment modeling. Files in the .sbi format can be converted
  to the .ssi format, which can store data about consecutive bases of the genome (segment). 
  This is useful when genomic context is important to some predictions, and could be used 
  as well for training models for CNV and structural rearrangements. This is work in progress.
* Draft SBI simulator (simulate-sbi.sh) to produce an SBI file corresponding to a VCF file, where
  counts are non-zero only for bases of a true genotype.
* Produce SBI files in a non-sorted format, suitable to convert to .ssi format (parallel-segments-sbi.sh)
* Draft SSI to SBI converter tool (sbi-to-ssi.sh). Takes a non-sorted .sbi file and produce a .ssi file where
  genotypes are represented base by base over a genomic segment.
  
### 1.3.2 (July 2017)

* parallel-gatk-realign.sh: add MD tags after realignment so that import to Goby or use with can proceed directly.
* parallel-gatk-realign.sh: optimizations to improve parallel processing speed.
* AddTrueGenotypeHelper (used by goby sequence_base_information output) fixed to correctly add indels.
* predict-genotype.sh: --no-cache option added to speed up one of predictions.
* new script: predict-genotypes-many.sh to process several sbi files in parallel 
  and call genotypes (to VCF format). This can be much faster than running goby to genotypes mode if you already have generated .sbi files.
* new genotype mapper: GenotypeMapperV37 has the best performance seen so far.
* new somatic mapper: SomaticFeatureMapper2, delegates to two copies of GenotypeMapperV37.
* write predicted somatic allele to VCF output: fixed org.campagnelab.dl.somatic.learning.architecture.graphs.SixDenseLayersNarrower2WithFrequencyAndBase
* TrainModel: Stop as soon as a NaN score is encountered.
* new script: generate-somatic-mutaed-set.sh, to generate somatic training set from two alignments.
* parallel-somatic-calls.sh: add focus varmap option to restrict output to specific sites in a varmap.
