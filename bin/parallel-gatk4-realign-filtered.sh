#!/usr/bin/env bash
if [ $# -lt 6 ]; then
 echo "usage: parallel-gatk4-realign-filtered.sh PATH_TO_GATK_LAUNCH 10g NUM_THREADS GENOME_FA BAM_INPUT BAM_OUTPUT [GATK_ARGS]"
 echo "usage: providing the an indexed bam (extension .bam.bai) next to the alignment will improve performance drastically."
 exit 1;
fi

GATK_LAUNCH=$1
MEMORY_PER_THREAD=$2
NUM_THREADS=$3
GENOME_FA=$4
BAM_INPUT=$5
BAM_OUTPUT=$6
GATK_ARGS=$7

DATE=`date +%Y-%m-%d`
echo ${DATE} >DATE.txt


echo "Using ${GATK_LAUNCH} to run GATK4."
echo "Using ${MEMORY_PER_THREAD} memory per thread."
echo "Using ${NUM_THREADS} number of threads."
echo "Using ${GENOME_FA} fasta file as genome. (index file .fa.fai required)."
echo "Using ${BAM_INPUT} as bam input."
echo "Appending GATK_ARGS: ${GATK_ARGS} to all GATK HaplotypeCaller calls"
echo "Writing bam output to ${BAM_OUTPUT}."

if [ ! -f ${BAM_INPUT}.bai ]; then
    samtools index ${BAM_INPUT}
fi

samtools idxstats ${BAM_INPUT} | cut -f 1 | head -n -1 |grep -v GL | sort -ur > refs.txt

# Disable read filters and over filters (needed for the model to see all features, bad and good):
NEEDED_OPTIONS="--max-reads-per-alignment-start 1000 --min-base-quality-score 0 --minimum-mapping-quality 0 \
--disable-read-filter MappingQualityReadFilter \
--disable-read-filter MappingQualityAvailableReadFilter \
--disable-read-filter NotSecondaryAlignmentReadFilter \
--disable-read-filter NotDuplicateReadFilter \
--disable-read-filter PassesVendorQualityCheckReadFilter \
--disable-read-filter NonZeroReferenceLengthAlignmentReadFilter \
--disable-read-filter GoodCigarReadFilter \
--disable-read-filter WellformedReadFilter"

--max-reads-per-alignment-start 1000 --min-base-quality-score 0 --minimum-mapping-quality 0 "
rm -rf calmd-and-convert-commands.txt
nLine=0
cat refs.txt | while read -r line
    do
       echo "${GATK_LAUNCH}  HaplotypeCaller --java-options \"-Xmx${MEMORY_PER_THREAD}\" \
            --spark-runner LOCAL -L ${line} --reference ${GENOME_FA} \
            --input ${BAM_INPUT} \
            ${GATK_ARGS} --output hc_variants_${nLine}.vcf  --bam-output realigned_slice_${nLine}.bam \
             ${NEEDED_OPTIONS} && \
         rm -f hc_variants_${nLine}.vcf* && \
         rm -f slice_${nLine}.bam && \
         rm -f slice_${nLine}.bam.bai && \
         samtools calmd -b -E realigned_slice_${nLine}.bam ${GENOME_FA} > realigned_slice_md_${nLine}.bam 2>&- && \
         samtools index  realigned_slice_md_${nLine}.bam && \
         rm -f realigned_slice_${nLine}.bam && \
         rm -f realigned_slice_${nLine}.bam.bai " >> calmd-and-convert-commands.txt
       nLine=$((nLine+1))
done

parallel --bar -j${NUM_THREADS} --progress --eta :::: calmd-and-convert-commands.txt
rm -fr ${BAM_OUTPUT}
samtools merge ${BAM_OUTPUT} realigned_slice_md_*.bam &&
samtools index ${BAM_OUTPUT}
rm realigned_slice_md_*.bam*