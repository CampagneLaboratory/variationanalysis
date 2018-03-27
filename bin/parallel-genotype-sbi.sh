#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

assertGobyInstalled
assertParallelInstalled


ALIGNMENTS="$*"
if [ -z "${OUTPUT_BASENAME+set}" ]; then

    if [ "$#" -eq 1 ]; then
       case ${ALIGNMENTS} in *.bam) OUTPUT_BASENAME=`basename ${ALIGNMENTS} .bam`;; esac
       case ${ALIGNMENTS} in *.sam) OUTPUT_BASENAME=`basename ${ALIGNMENTS} .sam`;; esac
       case ${ALIGNMENTS} in *.cram) OUTPUT_BASENAME=`basename ${ALIGNMENTS} .cram`;; esac
    else
        OUTPUT_BASENAME="out-concat"
    fi

    echo "OUTPUT_BASENAME set to ${OUTPUT_BASENAME}. Change the variable to influence the output basename."
fi

echo "Will write Goby alignment to ${OUTPUT_BASENAME}"
if [ -z "${SBI_GENOTYPE_VARMAP+set}" ]; then
    echo "Set SBI_GENOTYPE_VARMAP to the filename of a varmap to annotate true genotypes while you generate the .sbi file."
    VARMAP_OPTION=" "
else
   if [ -z "${REF_SAMPLING_RATE+set}" ]; then
        REF_SAMPLING_RATE="0.1"
        echo "Define REF_SAMPLING_RATE to change the default sampling rate. Set to ${REF_SAMPLING_RATE}"
    fi
    VARMAP_OPTION=" -x SequenceBaseInformationOutputFormat:random-seed=3764 -x SequenceBaseInformationOutputFormat:sampling-rate=${REF_SAMPLING_RATE} -x SequenceBaseInformationOutputFormat:true-genotype-map=${SBI_GENOTYPE_VARMAP} "
fi

if [ -z "${GOBY_NUM_SLICES+set}" ]; then
    GOBY_NUM_SLICES="50"
    echo "GOBY_NUM_SLICES set to ${GOBY_NUM_SLICES}. Change the variable to influence how many slices to use when processing an alignment."
fi
if [ -z "${INCLUDE_INDELS+set}" ]; then
    INCLUDE_INDELS="false"
    echo "INCLUDE_INDELS set to ${INCLUDE_INDELS}. Change the variable to influence whether indels are included in the SBI file."
fi
if [ -z "${SBI_GENOME+set}" ]; then
    SBI_GENOME="/data/genomes/Homo_sapiens.ucsc.hg19"
    echo "SBI_GENOME set to ${SBI_GENOME}. Change the variable to influence the genome used (must be indexed with goby build-sequence-cache)."
fi
if [ -z  "${SBI_NUM_THREADS+set}" ]; then
    SBI_NUM_THREADS="2"
    echo "SBI_NUM_THREADS set to ${SBI_NUM_THREADS}. Change the variable to influence the number of parallel jobs."
fi

if [ -z "${DO_CONCAT+set}" ]; then
    DO_CONCAT="true"
    echo "DO_CONCAT set to true."
fi

if [ -z "${DSV_OPTIONS+set}" ]; then
    DSV_OPTIONS="-n 0 -t 1"
fi
echo "DSV_OPTIONS set to ${DSV_OPTIONS}."

if [ -z "${REALIGN_AROUND_INDELS+set}" ]; then
    echo "Set REALIGN_AROUND_INDELS to true to enable realignment around indels."
    REALIGNMENT_OPTION=" "
else
   if [ "${REALIGN_AROUND_INDELS}" == "true" ]; then
       echo "REALIGN_AROUND_INDELS set to ${REALIGN_AROUND_INDELS}. ENABLED realignment around indels."
       REALIGNMENT_OPTION="--processor realign_near_indels"
   else
       echo "REALIGN_AROUND_INDELS set to ${REALIGN_AROUND_INDELS}. DISABLED realignment around indels."
       REALIGNMENT_OPTION=" "
   fi
fi

num_contigs=`goby 2g cfs ${ALIGNMENTS} --header-only|grep  'Number of target sequences'|head -1 | cut -d " " -f 6`
if [ ${GOBY_NUM_SLICES} -lt ${num_contigs} ]; then
    # make sure num slices is at least twice the number of chromosomes//contigs in the reference:
    GOBY_NUM_SLICES=`echo $(( num_contigs * 2 ))`
fi

goby 8g suggest-position-slices ${ALIGNMENTS} --modulo 1000 --number-of-slices ${GOBY_NUM_SLICES} --restrict-per-chromosome  -o slices.tsv
grep -v targetIdStart slices.tsv >slices


echo " discover-sequence-variants ${DSV_OPTIONS} --genome  ${SBI_GENOME} --format  SEQUENCE_BASE_INFORMATION  ${ALIGNMENTS} \
    --call-indels  ${INCLUDE_INDELS} ${REALIGNMENT_OPTION} \
    --max-coverage-per-site 1000 -x HTSJDKReaderImpl:force-sorted=true \
    -x SequenceBaseInformationOutputFormat:genomic-context-length=41 \
    ${VARMAP_OPTION} " >command.txt

# keep a log of the commands that were used to generate this dataset:
cp command.txt command-`date +%h_%d_%H_%M`.txt

cut -f3,6 slices  | awk 'BEGIN{count=1} {print "-s "$1" -e " $2" -o out-part-"(count++)}' >boundaries
cat boundaries
parallel --bar --eta -j${SBI_NUM_THREADS} --plus  --progress goby 10g  `cat command.txt`  :::: boundaries

# keep a log of the commands that were used to generate this dataset:
cp command.txt command-`date +%h_%d_%H_%M`.txt

# If the CHR variable is not set, use the default: "chr". This works when the genome has chromosomal contigs prefixed with "chr".
if [ -z "${CHR+set}" ]; then
    CHR="chr"
    echo "Setting CHR to 'chr'"
fi

if [ "${DO_CONCAT}" == "true" ]; then
    cat boundaries| grep -v -e "${CHR}19" -e "${CHR}20" -e "${CHR}21" -e "${CHR}22" -e "${CHR}X" -e "${CHR}Y" | cut -d " " ""-f 6 |awk '{print $1".sbi"}' >training-parts
    cat boundaries| grep -e "${CHR}19"  |cut -d " " ""-f 6 | awk '{print $1".sbi"}' >validation-parts
    cat boundaries| grep -v -e "${CHR}19" | grep -e "${CHR}20" -e "${CHR}21" -e "${CHR}22" -e "${CHR}X" -e "${CHR}Y" |cut -d " " ""-f 6 | awk '{print $1".sbi"}' >testing-parts

    if [ -s training-parts ]; then
        concat.sh ${memory_requirement} -f -i `cat training-parts`  -o ${OUTPUT_BASENAME}-pre-train
    fi
    if [ -s validation-parts ]; then
        concat.sh ${memory_requirement} -f -i `cat validation-parts` -o ${OUTPUT_BASENAME}-pre-validation
    fi
    if [ -s testing-parts ]; then
        concat.sh ${memory_requirement} -f -i `cat testing-parts` -o ${OUTPUT_BASENAME}-test
    fi
fi