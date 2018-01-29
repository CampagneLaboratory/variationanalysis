#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/common.sh
if [ "$#" -ne 2 ]; then
   echo "Argument missing. expected arguments memory_size goby_alignment goby_genome"
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
GENOME=$2
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

if [ -z  "${SBI_NUM_THREADS+set}" ]; then
    SBI_NUM_THREADS="2"
    echo "SBI_NUM_THREADS set to ${SBI_NUM_THREADS}. Change the variable to influence the number of parallel jobs."
fi

export SBI_GENOME=${GENOME}
rm -rf tmp
mkdir -p tmp

export SBI_GENOTYPE_VARMAP=" "
export REF_SAMPLING_RATE=1.0


set -x
# put all results under tmp:
export OUTPUT_BASENAME=tmp/${OUTPUT_PREFIX}

# restrict unlabeled sites to those that show some level of variations from the reference:
export DSV_OPTIONS="-n 5 -t 1"
unset SBI_GENOTYPE_VARMAP
DO_CONCAT="false"

parallel-genotype-sbi.sh 10g ${ALIGNMENTS} 2>&1 | tee parallel-genotype-sbi.log
dieIfError "Failed to generate .sbi file"
rm randomize-commands.txt

for file in out-part-*.sbi; do
        SBI_basename=`basename $file .sbi`
        if [ -f ${SBI_basename}.sbip ]; then
            echo "SBI basename: '$SBI_basename'"
            echo >> randomize-commands.txt "randomize.sh ${memory_requirement} -i ${file} -o random-${SBI_basename}.sbi  -b 100000  --random-seed 2378237"
        fi
done

parallel --bar --eta -j${SBI_NUM_THREADS} --plus  --progress :::: randomize-commands.txt

rm -rf tmp

concat.sh ${memory_requirement} -f -o ${OUTPUT_BASENAME}-unlabeled.sbi  -i  random-out-part-*.sbi



export DATASET="${OUTPUT_BASENAME}-"

# restrict unlabeled sites to those that show some level of variations from the reference:
export DSV_OPTIONS="-n 5 -t 1"
unset SBI_GENOTYPE_VARMAP
DO_CONCAT="false"

# parallel-genotype-sbi.sh 10g ${ALIGNMENTS} 2>&1 | tee parallel-genotype-sbi.log
# dieIfError "Failed to generate .sbi file"
rm randomize-commands.txt

for file in out-part-*.sbi; do
        SBI_basename=`basename $file .sbi`
        if [ -f ${SBI_basename}.sbip ]; then
            echo "SBI basename: '$SBI_basename'"
            echo >> randomize-commands.txt "randomize.sh ${memory_requirement} -i ${file} -o random-${SBI_basename}.sbi  -b 100000  --random-seed 2378237"
        fi
done

parallel --bar --eta -j${SBI_NUM_THREADS} --plus  --progress :::: randomize-commands.txt

rm -rf tmp

concat.sh ${memory_requirement} -f -o ${OUTPUT_BASENAME}-unlabeled  -i  random-out-part-*.sbi



export DATASET="${OUTPUT_BASENAME}-"
