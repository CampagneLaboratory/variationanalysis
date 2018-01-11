#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

assertGobyInstalled
assertParallelInstalled

if [ -e configure.sh ]; then
 echo "Loading configure.sh"
 source configure.sh
fi

ALIGNMENTS="$*"
if [ "$#" -eq 1 ]; then
   case ${ALIGNMENTS} in *.bam) OUTPUT_BASENAME=`basename ${ALIGNMENTS} .bam`;; esac
   case ${ALIGNMENTS} in *.sam) OUTPUT_BASENAME=`basename ${ALIGNMENTS} .sam`;; esac
   case ${ALIGNMENTS} in *.cram) OUTPUT_BASENAME=`basename ${ALIGNMENTS} .cram`;; esac
else
    OUTPUT_BASENAME="out-concat"
fi
echo "Will write Goby alignment to ${OUTPUT_BASENAME}"

if [ -z "${SBI_GENOME+set}" ]; then
    SBI_GENOME="/data/genomes/Homo_sapiens.ucsc.hg19"
    echo "SBI_GENOME set to ${SBI_GENOME}. Change the variable to influence the genome used (must be indexed with goby build-sequence-cache)."
fi
if [ -z "${FASTA_GENOME+set}" ]; then
    FASTA_GENOME="/data/genomes/Homo_sapiens.ucsc.hg19.fa"
    echo "FASTA_GENOME set to ${FASTA_GENOME}. Change the variable to influence the fasta used. Must come with a fasta index using same basename."
fi
if [ -z  "${SBI_NUM_THREADS+set}" ]; then
    SBI_NUM_THREADS="2"
    echo "SBI_NUM_THREADS set to ${SBI_NUM_THREADS}. Change the variable to influence the number of parallel jobs."
fi
echo "variables: ${SBI_GENOME} ${SBI_NUM_THREADS}"


set -x

samtools idxstats ${ALIGNMENTS} | cut -f 1 | head -n -1 > refs.txt
rm -rf calmd-and-convert-commands.txt

MD_LINES=`samtools view  ${ALIGNMENTS}  |head -1000|grep MD:Z:|wc -l`
if [ ${MD_LINES} -eq 1000 ]; then
    HAS_MD_TAGS="TRUE"

else
    HAS_MD_TAGS="FALSE"
fi

cat refs.txt | while read -r line
    do
       if [ ${HAS_MD_TAGS} == "TRUE" ]; then
           echo "\
             samtools view -u ${ALIGNMENTS} ${line} > md_${line}.bam 2> /dev/null && \
             samtools index md_${line}.bam 2>/dev/null &&\
             goby 8g concatenate-alignments --genome  ${SBI_GENOME}  md_${line}.bam  -o goby_slice_${line} &&\
             rm md_${line}.bam  &&\
             rm md_${line}.bam.bai \
           " >> calmd-and-convert-commands.txt

       else

           echo "\
           samtools calmd -E -u <(samtools view -b ${ALIGNMENTS} ${line}) ${FASTA_GENOME} > md_${line}.bam 2> /dev/null && \
             samtools index md_${line}.bam 2>/dev/null &&\
             goby 8g concatenate-alignments --genome  ${SBI_GENOME}  md_${line}.bam  -o goby_slice_${line} &&\
             rm md_${line}.bam  &&\
             rm md_${line}.bam.bai \
           " >> calmd-and-convert-commands.txt
        fi
done

parallel  --results error-folder.log --bar -j${SBI_NUM_THREADS} --eta :::: calmd-and-convert-commands.txt

goby ${memory_requirement} concatenate-alignments goby_slice_*.entries -o ${OUTPUT_BASENAME} | tee concat.log && rm goby_slice_*