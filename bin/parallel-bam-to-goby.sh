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
    SBI_GENOME="/data/genomes/Homo_sapiens.ucsc.hg19.fa"
    echo "SBI_GENOME set to ${SBI_GENOME}. Change the variable to influence the fasta used."
fi
if [ -z  "${SBI_NUM_THREADS+set}" ]; then
    SBI_NUM_THREADS="2"
    echo "SBI_NUM_THREADS set to ${SBI_NUM_THREADS}. Change the variable to influence the number of parallel jobs."
fi
echo "variables: ${SBI_GENOME} ${SBI_NUM_THREADS}"

goby ${memory_requirement} suggest-position-slices ${ALIGNMENTS} --number-of-slices 60-o slices.tsv
grep -v targetIdStart slices.tsv >slices

nLine=0
tail -n +2 file.txt | while read line
    do
       sRef=`cut -f1 line`
       sPos=`cut -f2 line`
       ePos=`cut -f5 line`
       samtools view -b ${ALIGNMENTS} ${sReg}:${sPos}-${ePos} > slice_${nLine}.bam
       echo "samtools calmd $slice_{nLine}.bam ${FASTA_GENOME} > md_slice_${nLine}.bam ;\
        goby ${memory_requirement} concatenate-alignments --genome  ${SBI_GENOME}  md_slice_${nLine}.bam -o goby_slice_${nLine}"\
         >> calmd-and-convert-commands.txt
       var=$((var+1))
       fi
    done
done

for
parallel -j${SBI_NUM_THREADS} -a --progress calmd-and-convert-commands.txt echo

goby ${memory_requirement} concatenate-alignments goby_slice_*.entries -o ${OUTPUT_BASENAME}
