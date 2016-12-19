Starting with release 1.2, we are offering a preview of this module.
See the forum at the bottom for questions and our call for collaborators.

The genotype module provides tools to call genotypes with high-throughput
sequencing data are available as well as the ground-truth necessary to
train the models.

Genotype offers the following tools:

 - add-true-genotypes.sh. This tool reads a .sbi file produced with
 Goby (see [tutorial](../GENOTYPE-TUTORIAL.md)) and adds ground-truth
 genotypes. It generates an annotated .sbi file that can be used for
 training or evaluating models. Ground-truth must be provided in a binary
 map file which can also be produced by Goby, as shown in the tutorial.

- parallel-genotype-sbi.sh. This tool uses GNU parallel to speed up
the creation of an .sbi file for an alignment (in Goby, BAM or CRAM format).
It splits up
 the input alignment(s) into genomic slices and converts each slice
 independently. Slices .sbi are combined with concat.sh. The tool accepts
 as argument the input alignments. It assumes alignment against hg19, but
 you can change this by defining the SBI_GENOME environment variable
 to point to the basename of the Goby/samtools indexed genome.

 export SBI_GENOME=path/to/genome.
 The variable SBI_NUM_THREADS controls how many parallel threads will be
 used. Reads the script in distribution/bin/ to see the default values.

- train-genotype.sh. This tool trains models for genotypes. See the
tutorial and documentation for its various options.

- predict-genotypes.sh. This tool accepts a trained model and an input
.sbi file (annotated or not) and predicts the genotype of each site in the
.sbi file. If ground-truth annotations  are provided, the output indicates
if predictions are correct or wrong and provide options to filter on correctness
and probability. In combination with show-genotypes.sh, this is a useful
tool for error analysis.

- show-genotypes.sh. This tool can show more details about genotype
predictions done with predict-genotypes.sh.


For more information, see the [project forum](https://groups.google.com/forum/#!forum/variationanalysis).