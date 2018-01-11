#!/usr/bin/env bash
if [ ! "$#" -eq 3 ]; then
    echo "Usage extract-chromosomes.sh NUM_THREADS BAM_FILE CHROMOSOME_LIST_SPACE_SEPARATED"
    exit 1
fi

NUM_THREADS=$1
BAM_FILE=$2
CHROMOSOME_LIST=$3
BASENAME=`basename ${BAM_FILE} .bam`

for CHROMOSOME in ${CHROMOSOME_LIST}; do
 echo "Extracting $CHROMOSOME"
 sambamba slice -o "${BASENAME}-${CHROMOSOME}.bam" ${BAM_FILE} ${CHROMOSOME}
 sambamba index -t ${NUM_THREADS} "${BASENAME}-${CHROMOSOME}.bam"
done

rm -fr ${BASENAME}-"subset.bam"

if [ `echo ${CHROMOSOME_LIST} |wc -w ` -gt 1 ]; then
    echo "Merging"
    sambamba merge -t ${NUM_THREADS} ${BASENAME}-"subset.bam" ${BASENAME}-*.bam
else
    mv ${BASENAME}-${CHROMOSOME_LIST}.bam ${BASENAME}-"subset.bam"
fi

echo "Indexing"
sambamba index -t ${NUM_THREADS} ${BASENAME}-"subset.bam"

echo "Done"
