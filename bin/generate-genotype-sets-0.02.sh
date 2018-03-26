#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/common.sh
if [ "$#" -ne 4 ]; then
   echo "Argument missing. expected arguments memory_size goby_alignment vcf goby_genome"
   exit 1;
fi

if [ -e configure.sh ]; then
 echo "Loading configure.sh"
 source configure.sh
fi

memory_requirement=$1
#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh
ALIGNMENTS=$1
VCF=$2
GENOME=$3
DATE=`date +%Y-%m-%d`
echo ${DATE} >DATE.txt

case ${ALIGNMENTS} in *.bam) OUTPUT_PREFIX=`basename ${ALIGNMENTS} .bam`;; esac
case ${ALIGNMENTS} in *.sam) OUTPUT_PREFIX=`basename ${ALIGNMENTS} .sam`;; esac
case ${ALIGNMENTS} in *.cram) OUTPUT_PREFIX=`basename ${ALIGNMENTS} .cram`;; esac
case ${ALIGNMENTS} in *.entries) OUTPUT_PREFIX=`basename ${ALIGNMENTS} .entries`;; esac
OUTPUT_PREFIX="${OUTPUT_PREFIX}-${DATE}"
echo "Will write results to ${OUTPUT_PREFIX}"

if [ -z "${DELETE_TMP}" ]; then
    DELETE_TMP="false"
    echo "DELETE_TMP set to ${DELETE_TMP}. Change the variable with export to clear the working directory."
fi

if [ -z "${SBI_SPLIT_OVERRIDE_DESTINATION+set}" ]; then
   export SBI_SPLIT_OVERRIDE_DESTINATION_OPTION=""
   echo "SBI_SPLIT_OVERRIDE_DESTINATION not set. Change the variable to chr21,chr22 to include only chr21 and chr22 in the test set."
else
    export SBI_SPLIT_OVERRIDE_DESTINATION_OPTION=" --destination-override test:${SBI_SPLIT_OVERRIDE_DESTINATION} "
    echo "Using SBI_SPLIT_OVERRIDE_DESTINATION=${SBI_SPLIT_OVERRIDE_DESTINATION} to put chromosomes ${SBI_SPLIT_OVERRIDE_DESTINATION} into test set."
fi

if [-z ${VARMAP_CHR_PREFIX} ]; then
    VARMAP_CHR_PREFIX_OPTION="--chromosome-prefix ${VARMAP_CHR_PREFIX}"
else
    VARMAP_CHR_PREFIX_OPTION=" "
fi

export SBI_GENOME=${GENOME}
rm -rf tmp
mkdir -p tmp
if [ -z "${SBI_GENOTYPE_VARMAP+set}" ]; then
  echo "Recalculating varmap. Set SBI_GENOTYPE_VARMAP to reuse a previous map."
  goby ${memory_requirement} vcf-to-genotype-map ${VCF} --genome ${SBI_GENOME} \
        -o ${OUTPUT_PREFIX}.varmap ${VARMAP_CHR_PREFIX_OPTION}
  dieIfError "Failed to produce varmap"
  echo "export SBI_GENOTYPE_VARMAP=${OUTPUT_PREFIX}.varmap" >> configure.sh
  . configure.sh
fi


if [ -z "${REF_SAMPLING_RATE+set}" ]; then
    # We keep only 10% of reference matching sites as we add the true genotypes:
    REF_SAMPLING_RATE="0.02"
    echo "REF_SAMPLING_RATE set to ${REF_SAMPLING_RATE}. Change the variable to influence what percent of reference sites are sampled."
else
 echo "REF_SAMPLING_RATE set to ${REF_SAMPLING_RATE}."
fi

set -x
# put all results under tmp:
export OUTPUT_BASENAME=tmp/${OUTPUT_PREFIX}
export DO_CONCAT="true"

parallel-genotype-sbi.sh 10g ${ALIGNMENTS} 2>&1 | tee parallel-genotype-sbi.log
dieIfError "Failed to generate .sbi file"

export OUTPUT_BASENAME=${OUTPUT_PREFIX}
if [ -s tmp/${OUTPUT_BASENAME}-pre-train.sbip ]; then
    # randomize and put final datasets in current directory:
    randomize.sh ${memory_requirement} -i tmp/${OUTPUT_BASENAME}-pre-train.sbi -o ${OUTPUT_BASENAME}-train.sbi  \
       -b 100000 -c 100   --random-seed 2378237 |tee randomize.log
    dieIfError "Failed to randomize training set."
fi
if [ -s tmp/${OUTPUT_BASENAME}-pre-validation.sbip ]; then
    randomize.sh ${memory_requirement} -i tmp/${OUTPUT_BASENAME}-pre-validation.sbi -o ${OUTPUT_BASENAME}-validation.sbi  \
     -b 100000 -c 100   --random-seed 2378237 |tee randomize.log
    dieIfError "Failed to randomize validation set"
fi
if [ -s tmp/${OUTPUT_BASENAME}-test.sbip ]; then
    mv tmp/${OUTPUT_BASENAME}-test.sbi* ./
    dieIfError "Failed to copy test set"
fi
# no need to randomize the test set.

if [ ${DELETE_TMP} = "true" ]; then
   rm tmp/${OUTPUT_PREFIX}-pre-training.sbi
   rm tmp/${OUTPUT_PREFIX}-pre-training.sbi
   rm -rf tmp
fi

export DATASET="${OUTPUT_BASENAME}-"
