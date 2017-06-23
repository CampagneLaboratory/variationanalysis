GATK_JAR=$1
MEMORY_PER_THREAD=$2
NUM_THREADS=$3
GENOME_FA=$4
KNOWN_INDELS_VCF=$5
BAM_INPUT=$6
BAM_OUTPUT=$7

DATE=`date +%Y-%m-%d`
echo ${DATE} >DATE.txt


echo "Using ${GATK_JAR} as gatk jar."
echo "Using ${MEMORY_PER_THREAD} memory per thread."
echo "Using ${NUM_THREADS} number of threads."
echo "Using ${GENOME_FA} fasta file as genome. (index file .fa.fai required)."
echo "Using ${KNOWN_INDELS_VCF} vcf file with known indel set."
echo "Using ${BAM_INPUT} as bam input."
echo "Writing bam output to ${BAM_OUTPUT}."

if [ ! -f ${BAM_INPUT}.bai ]; then
    samtools index ${BAM_INPUT}
fi

samtools idxstats ${BAM_INPUT} | cut -f 1 | head -n -1 > refs.txt

rm -rf calmd-and-convert-commands.txt
nLine=0
cat refs.txt | while read -r line
    do
       echo "samtools view -u ${BAM_INPUT} ${line} > slice_${nLine}.bam && \
         samtools index slice_${nLine}.bam && \
         java -jar -Xmx${MEMORY_PER_THREAD} ${GATK_JAR} -R ${GENOME_FA} -ip 50 -T RealignerTargetCreator \
                -known ${KNOWN_INDELS_VCF} -I slice_${nLine}.bam  -o slice_${nLine}_realignment_targets.interval_list \
                -mismatch 0.0 && \
         java -jar -Xmx${MEMORY_PER_THREAD} ${GATK_JAR} -R ${GENOME_FA} -ip 50 -T IndelRealigner \
         -known ${KNOWN_INDELS_VCF} -I slice_${nLine}.bam -targetIntervals slice_${nLine}_realignment_targets.interval_list \
                -o realigned_slice_${nLine}.bam && \
         rm  slice_${nLine}_realignment_targets.interval_list && \
         rm slice_${nLine}.bam && \
         rm slice_${nLine}.bam.bai \
         samtools calmd -E -u realigned_slice_${nLine}.bam ${GENOME_FA} > realigned_slice_md_${nLine}.bam 2>&- && \
         samtools index realigned_slice_md_${nLine}.bam  &&
         rm -f realigned_slice_${nLine}.bam && \
         rm -f realigned_slice_${nLine}.bam.bai
       " >> calmd-and-convert-commands.txt
       nLine=$((nLine+1))
done

parallel --bar -j${NUM_THREADS} --eta :::: calmd-and-convert-commands.txt

samtools merge -f ${BAM_OUTPUT} realigned_slice_md_*.bam &&

rm realigned_md_slice_*