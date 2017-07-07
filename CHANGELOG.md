##CHANGE LOG

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
