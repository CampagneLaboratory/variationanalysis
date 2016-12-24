if [ "$#" -ne 4 ]; then
   echo "Argument missing. expected arguments memory_size goby_alignment vcf goby_genome"
   exit 1;
fi

memory_requirement=$1
#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh
ALIGNMENTS=$1
VCF=$2
GENOME=$3


if [ -z "${DELETE_TMP}" ]; then
    DELETE_TMP="false"
    echo "DELETE_TMP set to ${DELETE_TMP}. Change the variable with export to clear the working directory."
fi
rm -rf tmp
mkdir tmp

goby ${memory_requirement} vcf-to-genotype-map ${VCF} \
  -o tmp/variants.varmap

export SBI_GENOME=${GENOME}
export OUTPUT_BASENAME=tmp/genotype_full

parallel-genotype-sbi.sh ${memory_requirement} ${ALIGNMENTS}


add-true-genotypes.sh ${memory_requirement} -m tmp/variants.varmap \
  -i tmp/genotype_full.sbi \
  -o tmp/genotype_full_called.sbi \
  --genome ${GENOME}

randomize.sh ${memory_requirement} -i tmp/genotype_full_called.sbi \
  -o tmp/genotype_full_called_randomized.sbi

split.sh ${memory_requirement} -i tmp/genotype_full_called_randomized.sbi \
  -f 0.8 -f 0.1 -f 0.1 \
  -o genotype_noIndel_ \
   -s train -s test -s validation

# subset the validation sample, throwing out many reference matching sites (to speed
# up performance evaluation for early stopping):
mv genotype_noIndel_validation.sbi genotype_noIndel_validation-all.sbi
mv genotype_noIndel_validation.sbip genotype_noIndel_validation-all.sbip

add-true-genotypes.sh ${memory_requirement} -m tmp/variants.varmap \
  -i genotype_noIndel_validation-all.sbi \
  -o genotype_noIndel_validation \
  --genome ${GENOME} --ref-sampling-rate 0.01

if [ ${DELETE_TMP} = "true" ]; then
   rm -rf tmp
fi

export DATASET=genotype_noIndel_