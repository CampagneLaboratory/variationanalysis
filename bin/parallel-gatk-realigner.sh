GATK_JAR=$0
MEMORY_PER_THREAD=$1
NUM_THREADS=$2
GENOME_FA=$3
BAM_INPUT=$4
BAM_OUTPUT=$5

DATE=`date +%Y-%m-%d`
echo ${DATE} >DATE.txt


echo "Using ${GATK_JAR} as gatk jar."
echo "Using ${MEMORY_PER_THREAD} memory per thread."
echo "Using ${NUM_THREADS} number of threads."
echo "Using ${GENOME_FA} fasta file as genome. (index file .fa.fai required)."
echo "Using ${BAM_INPUT} as bam input."
echo "Writing bam output to ${BAM_OUTPUT}."

echo "variables: ${SBI_GENOME} ${SBI_NUM_THREADS}"



if [! -f ${BAM_INPUT}.bai]; then
    samtools index ${BAM_INPUT}
if

goby 8g suggest-position-slices ${ALIGNMENTS} --number-of-slices 60 -o slices.tsv --restrict-per-chromosome

rm -rf calmd-and-convert-commands.txt
nLine=0
tail -n +2 slices.tsv | while read -r line
    do
       sRef=`echo $line | cut -f1 -d ' '`
       sPos=`echo $line | cut -f2 -d ' '`
       ePos=`echo $line | cut -f5 -d ' '`
       echo "samtools view -u ${ALIGNMENTS} ${sRef}:${sPos}-${ePos} > slice_${nLine}.bam ;\
         samtools calmd -E -u slice_${nLine}.bam ${FASTA_GENOME} > md_slice_${nLine}.bam ;\
         samtools index md_slice_${nLine}.bam &&\
         rm slice_${nLine}.bam ;\
         java -Xmx${MEMORY_PER_THREAD} -jar ${GATK_JAR} -T HaplotypeCaller -R ${GENOME_FA} -I md_slice_${nLine}.bam -o hc_variants.vcf  -bamout realigned_md_slice_${nLine}.bam &&\
         rm md_slice_${nLine}.bam &&\
         rm md_slice_${nLine}.bam.bai \
       " >> calmd-and-convert-commands.txt
       nLine=$((nLine+1))
done

parallel --bar -j${NUM_THREADS} --eta :::: calmd-and-convert-commands.txt

samtools merge ${BAM_OUTPUT} realigned_md_slice_*.entries &&

rm goby_slice_*