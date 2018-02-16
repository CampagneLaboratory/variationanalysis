#!/usr/bin/env bash

NUM_ARGS="$#"
NUM_ARGS_EXPECTED="${NUM_ARGS}"
. `dirname "${BASH_SOURCE[0]}"`/common.sh
if [ -z "${VCF_OUTPUT+set}" ] || [ -z "${BED_OBSERVED_REGIONS_OUTPUT+set}" ]; then
  NUM_ARGS_EXPECTED=5
fi

if [ "${NUM_ARGS}" == 5 ]; then
  unset VCF_OUTPUT
  unset BED_OBSERVED_REGIONS_OUTPUT
fi

if [ ! "${NUM_ARGS}" == "${NUM_ARGS_EXPECTED}" ]; then
   echo "Usage: evaluate-genotypes-vec.sh model-directory [checkpoint-key model-prefix test-sbi dataset]."
   echo "The env variables GOLD_STANDARD_VCF_SNP_GZ GOLD_STANDARD_VCF_INDEL_GZ and GOLD_STANDARD_CONFIDENT_REGIONS_BED_GZ can be used to change the VCFs and confident region bed."
   echo "The first run downloads these files from the Genome in a Bottle for sample NA12878 when the variables are not defined."
   echo "The first argument should point to the directory of the trained PyTorch model being used with a config.properties within."
   echo "The 2nd-5th arguments are optional. The 2nd argument should be the checkpoint key used and the 3rd should be the model prefix (usually best or latest)."
   echo "The 4th argument should be a path to the test SBI file, and the 5th should be the name of the dataset used for the sbi."
   echo "The SBI file should be named according to the convention <basename>-<dataset>.sbi, with a corresponding <basename>-<dataset>.vec file in the same directory."
   echo "You can bypass the predict phase by defining the variables VCF_OUTPUT and BED_OBSERVED_REGIONS_OUTPUT to point to the output of predict"
   exit 1;
fi
MODEL_DIR=$1

if [ -e configure.sh ]; then
 echo "Loading configure.sh"
 source configure.sh
fi

if [ -z "${RTG_TEMPLATE+set}" ]; then
  RTG_TEMPLATE="hg19.sdf"
  echo "RTG_TEMPLATE not set, using default=${RTG_TEMPLATE}"
fi

if [ ! -e "${RTG_TEMPLATE}" ]; then
 echo "You must install an rtg template, or build one (rtg format file.fa -o template.sdf) in the current directory. See rtg downloads at http://www.realtimegenomics.com/news/pre-formatted-reference-datasets/"
 exit 10;
fi

function assertRTGInstalled {
    echo no | rtg version >/dev/null 2>&1 || { echo >&2 "This script requires rtg but it's not installed. Aborting. Install rtg (see https://github.com/genome-in-a-bottle/giab_FAQ) and add the rtg executable to your path, then try again."; exit 1; }
}
assertRTGInstalled

function assertVcfToolsInstalled {
   echo done| vcf-sort >/dev/null 2>&1 || { echo >&2 "This script requires vcf-sort but it's not installed. Aborting. Install vcf-tools (see https://sourceforge.net/projects/vcftools/files/) and add the vcf-sort executable to your path, then try again."; exit 1; }
}

assertVcfToolsInstalled

function assertBedtoolsInstalled {
   bedtools --version  >/dev/null 2>&1 || { echo >&2 "This script requires bedtools but it's not installed. Aborting. Install bedtools (see http://bedtools.readthedocs.io/en/latest/content/installation.html) and add the bedtools executable to your path, then try again."; exit 1; }
}

assertBedtoolsInstalled

function assertNumpyInstalled {
    python -c "import numpy" >/dev/null 2>&1 || { echo >&2 "This script requires numpy to be installed. Aborting. Check to make sure numpy (http://www.numpy.org) is installed in your current environment or venv."; exit 1; }
}

assertNumpyInstalled

function assertTorchInstalled {
    python -c "import torch" >/dev/null 2>&1 || { echo >&2 "This script requires pytorch to be installed. Aborting. Check to make sure pytorch (http://pytorch.org) is installed in your current environment or venv."; exit 1; }
}

assertTorchInstalled

function assertTorchnetInstalled {
    python -c "import torchnet" >/dev/null 2>&1 || { echo >&2 "This script requires torchnet to be installed. Aborting. Check to make sure torchnet (https://github.com/pytorch/tnt) is installed in your current environment or venv."; exit 1; }
}

assertTorchnetInstalled

function assertPredictDatasetInstalled {
    predict-dataset.sh --help >/dev/null 2>&1 || { echo >&2 "This script requires predict-dataset.sh from GenotypeTensors to be in your path. Aborting. Check to make sure it is in your path before running, then try again."; exit 1; }
}

assertPredictDatasetInstalled

if [ -z "${MINI_BATCH_SIZE+set}" ]; then
    MINI_BATCH_SIZE="2048"
    echo "MINI_BATCH_SIZE set to ${MINI_BATCH_SIZE}. Change the variable to switch the mini-batch-size."
fi

set -x
if [ -z "${GOLD_STANDARD_VCF_SNP_GZ+set}" ] || [ -z "${GOLD_STANDARD_VCF_INDEL_GZ+set}" ]; then
    if [ -z "${GOLD_STANDARD_VCF_GZ+set}" ]; then
        echo "Downloading Gold standard Genome in a Bottle VCF. Define GOLD_STANDARD_VCF_GZ to use an alternate Gold Standard."
        rm -fr HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf_phased.vcf.gz.*
        wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.1/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf_phased.vcf.gz
        mv HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf_phased.vcf.gz GIAB-NA12878-confident.vcf.gz
        # add "chr prefix:"
        gzip -c -d GIAB-NA12878-confident.vcf.gz|awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' >GIAB-NA12878-confident-chr.vcf
        bgzip -f GIAB-NA12878-confident-chr.vcf
        tabix -f GIAB-NA12878-confident-chr.vcf.gz
        echo 'export GOLD_STANDARD_VCF_GZ="GIAB-NA12878-confident-chr.vcf.gz"' >>configure.sh
        export GOLD_STANDARD_VCF_GZ="GIAB-NA12878-confident-chr.vcf.gz"
    else
      echo "Formatting GOLD_STANDARD_VCF_GZ VCF for SNPs and indels"
    fi

    # remove non-SNPs:
    gzip -c -d  ${GOLD_STANDARD_VCF_GZ} |awk '{if($0 !~ /^#/) { if (length($4)==1 && length($5)==1) print $0;}  else {print $0}}' >GOLD-confident-chr-snps.vcf
    bgzip -f GOLD-confident-chr-snps.vcf
    tabix -f GOLD-confident-chr-snps.vcf.gz
    export GOLD_STANDARD_VCF_SNP_GZ="GOLD-confident-chr-snps.vcf.gz"

    # keep only indels:
    gzip -c -d  ${GOLD_STANDARD_VCF_GZ} |awk '{if($0 !~ /^#/) { if (length($4)!=1 || length($5)!=1) print $0;}  else {print $0}}' >GOLD-confident-chr-indels.vcf
    bgzip -f GOLD-confident-chr-indels.vcf
    tabix -f GOLD-confident-chr-indels.vcf.gz
    export GOLD_STANDARD_VCF_INDEL_GZ="GOLD-confident-chr-indels.vcf.gz"
    echo 'export GOLD_STANDARD_VCF_SNP_GZ="GOLD-confident-chr-snps.vcf.gz"' >>configure.sh
    echo 'export GOLD_STANDARD_VCF_INDEL_GZ="GOLD-confident-chr-indels.vcf.gz"' >>configure.sh
    echo "Gold standard VCF downloaded for NA12878 (SNPs) and named in configure.sh. Edit GOLD_STANDARD_VCF_SNP_GZ/GOLD_STANDARD_VCF_INDEL_GZ to switch to a different gold-standard validation VCF."
fi


if [ -z "${GOLD_STANDARD_CONFIDENT_REGIONS_BED_GZ+set}" ]; then
    echo "Downloading Gold standard Genome in a Bottle Confident Regions (bed)"
    rm -fr HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf.bed*
    wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.1/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf.bed
    mv HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf.bed GIAB-NA12878-confident-regions.bed
    # add "chr prefix:"
    cat  GIAB-NA12878-confident-regions.bed |awk '{print "chr"$1"\t"$2"\t"$3}' >GIAB-NA12878-confident-regions-chr.bed
    bgzip -f GIAB-NA12878-confident-regions-chr.bed
    tabix -f GIAB-NA12878-confident-regions-chr.bed.gz
    rm GIAB-NA12878-confident-regions.bed

    GOLD_STANDARD_CONFIDENT_REGIONS_BED_GZ="GIAB-NA12878-confident-regions-chr.bed.gz"
    echo 'export GOLD_STANDARD_CONFIDENT_REGIONS_BED_GZ="GIAB-NA12878-confident-regions-chr.bed.gz"' >>configure.sh
    echo "Gold standard confident regions downloaded for NA12878  and named in configure.sh. Edit GOLD_STANDARD_CONFIDENT_REGIONS_BED_GZ to switch to a different gold-standard confident region bed file."
fi

if [ "${NUM_ARGS}" == 5 ]; then
    CHECKPOINT_KEY=$2
    MODEL_PREFIX=$3
    DATASET_SBI=$4
    if [ ! -e "${DATASET_SBI}" ]; then
        echo "The SBI file was not found: ${DATASET_SBI}  "
        exit 1;
    fi
    DATASET_NAME=$5
    DATASET_BASENAME="`dirname ${DATASET_SBI}`/`basename ${DATASET_SBI} "-${DATASET_NAME}.sbi"`"
    DATASET_VEC="${DATASET_BASENAME}-${DATASET_NAME}.vec"
    if [ ! -e "${DATASET_VEC}" ]; then
        echo "The vector file was not found: ${DATASET_VEC}"
        exit 1;
    fi
fi


if [ -z "${MAX_RECORDS+set}" ]; then
  export TRAIN_MAX_RECORDS=""
  export PREDICT_MAX_RECORDS=""
  echo "MAX_RECORDS not set. Use it to limit how many records are used to training, validation and testing. Set to 10,000 or 100,000 for quick iterations."
else
export TRAIN_MAX_RECORDS="-n ${MAX_RECORDS} -x ${MAX_RECORDS}"
export PREDICT_MAX_RECORDS="-n ${MAX_RECORDS}"
fi
PREDICT_OPTIONS=""

if [ -z "${VCF_OUTPUT+set}" ] || [ -z "${BED_OBSERVED_REGIONS_OUTPUT+set}" ]; then
    echo "VCF_OUTPUT or BED_OBSERVED_REGIONS_OUTPUT are not defined. Running predict for ${DATASET_SBI}."
    MODEL_TIME=`basename ${MODEL_DIR}`
    echo "Running predict-dataset to create output vector file from trained PyTorch model..."
    predict-dataset.sh --model-path ${MODEL_DIR} --model-label ${MODEL_PREFIX} \
        --problem "genotyping:${DATASET_BASENAME}" \
        --dataset validation --mini-batch-size ${MINI_BATCH_SIZE} \
        --output "${DATASET_BASENAME}_predicted" --checkpoint-key ${CHECKPOINT_KEY}
    dieIfError "Failed to predict dataset with model ${MODEL_DIR}/."

    echo "Running predict-genotypes to create VCF and observed region bed.."
    predict-genotypes.sh 20g -m ${MODEL_DIR} -l ${MODEL_PREFIX} -f -i ${DATASET_SBI} \
        --format VCF --mini-batch-size ${MINI_BATCH_SIZE} --vec-path "${DATASET_BASENAME}_predicted.vec" \
         ${PREDICT_OPTIONS}  --no-cache ${PREDICT_MAX_RECORDS}
    dieIfError "Failed to create vcf from ${MODEL_DIR}/."
    echo "Evaluation with rtg vcfeval starting.."

    export VCF_OUTPUT=`ls -1tr ${MODEL_TIME}-${MODEL_PREFIX}-*.vcf|tail -1`
    export BED_OBSERVED_REGIONS_OUTPUT=`ls -1tr ${MODEL_TIME}-${MODEL_PREFIX}-*-observed-regions.bed |tail -1`
else
    echo "Evaluating with VCF_OUTPUT=${VCF_OUTPUT} and BED_OBSERVED_REGIONS_OUTPUT=${BED_OBSERVED_REGIONS_OUTPUT}"
fi

export VCF_OUTPUT_SORTED=`basename ${VCF_OUTPUT} .vcf`-sorted.vcf

if [ ! -e "${VCF_OUTPUT_SORTED}.gz" ]; then


    cat ${VCF_OUTPUT} | vcf-sort > ${VCF_OUTPUT_SORTED}
    dieIfError "Unable to sort prediction VCF."

    bgzip -f ${VCF_OUTPUT_SORTED}
    tabix ${VCF_OUTPUT_SORTED}.gz
fi

if [ ! -e "${BED_OBSERVED_REGIONS_OUTPUT}-sorted.bed.gz" ]; then

    # note -V option on chromosome key below is necessary on Centos 7, with sort version 8.22,
    # see https://github.com/chapmanb/bcbio-nextgen/issues/624
    sort -k1,1V -k2,2n ${BED_OBSERVED_REGIONS_OUTPUT} | bedtools merge > ${BED_OBSERVED_REGIONS_OUTPUT}-sorted.bed
    bgzip -f ${BED_OBSERVED_REGIONS_OUTPUT}-sorted.bed
    tabix ${BED_OBSERVED_REGIONS_OUTPUT}-sorted.bed.gz
fi


RTG_OUTPUT_FOLDER=output-${RANDOM}
gzip -c -d ${VCF_OUTPUT_SORTED}.gz |awk '{if($0 !~ /^#/) { if (length($4)==1 && length($5)==1) print $0;}  else {print $0}}'  >${VCF_OUTPUT_SORTED}-snps.vcf
bgzip -f ${VCF_OUTPUT_SORTED}-snps.vcf
tabix -f ${VCF_OUTPUT_SORTED}-snps.vcf.gz

rtg vcfeval --baseline=${GOLD_STANDARD_VCF_SNP_GZ}  \
        -c ${VCF_OUTPUT_SORTED}-snps.vcf.gz -o ${RTG_OUTPUT_FOLDER}/snp --template=${RTG_TEMPLATE}  \
            --evaluation-regions=${GOLD_STANDARD_CONFIDENT_REGIONS_BED_GZ} \
            --bed-regions=${BED_OBSERVED_REGIONS_OUTPUT}-sorted.bed.gz \
            --vcf-score-field=P  --sort-order=descending
dieIfError "Failed to run rtg vcfeval for SNPs."

cp ${VCF_OUTPUT_SORTED}-snps.vcf.gz  ${RTG_OUTPUT_FOLDER}/snp/

gzip -c -d ${VCF_OUTPUT_SORTED}.gz |awk '{if($0 !~ /^#/) { if (length($4)!=1 || length($5)!=1) print $0;}  else {print $0}}'  >${VCF_OUTPUT_SORTED}-indels.vcf
bgzip -f ${VCF_OUTPUT_SORTED}-indels.vcf
tabix -f ${VCF_OUTPUT_SORTED}-indels.vcf.gz

rtg vcfeval --baseline=${GOLD_STANDARD_VCF_INDEL_GZ}  \
        -c ${VCF_OUTPUT_SORTED}-indels.vcf.gz -o ${RTG_OUTPUT_FOLDER}/indel --template=${RTG_TEMPLATE}  \
            --evaluation-regions=${GOLD_STANDARD_CONFIDENT_REGIONS_BED_GZ} \
            --bed-regions=${BED_OBSERVED_REGIONS_OUTPUT}-sorted.bed.gz \
            --vcf-score-field=P  --sort-order=descending
dieIfError "Failed to run rtg vcfeval."

cp ${VCF_OUTPUT_SORTED}-indels.vcf.gz  ${RTG_OUTPUT_FOLDER}/indel/

MODEL_TIME=`basename ${MODEL_DIR}`
grep ${MODEL_TIME} predict-statistics.tsv   >${RTG_OUTPUT_FOLDER}/predict-statistics.tsv


cp ${MODEL_DIR}/config.properties ${RTG_OUTPUT_FOLDER}
cp ${BED_OBSERVED_REGIONS_OUTPUT}-sorted.bed.gz ${RTG_OUTPUT_FOLDER}
echo "See rtg vcfeval detailed output in ${RTG_OUTPUT_FOLDER}"

RTG_ROCPLOT_OPTIONS="--scores"
rtg rocplot ${RTG_OUTPUT_FOLDER}/snp/snp_roc.tsv.gz --svg ${RTG_OUTPUT_FOLDER}/snp/SNP-ROC.svg ${RTG_ROCPLOT_OPTIONS} --title="SNPs, model ${MODEL_TIME}"
dieIfError "Unable to generate SNP ROC plot."

rtg rocplot ${RTG_OUTPUT_FOLDER}/snp/snp_roc.tsv.gz -P --svg ${RTG_OUTPUT_FOLDER}/snp/SNP-PrecisionRecall.svg ${RTG_ROCPLOT_OPTIONS} --title="SNPs, model ${MODEL_TIME}"
dieIfError "Unable to generate SNP Precision Recall plot."

RTG_ROCPLOT_OPTIONS="--scores"
rtg rocplot ${RTG_OUTPUT_FOLDER}/indel/non_snp_roc.tsv.gz -P --svg ${RTG_OUTPUT_FOLDER}/indel/INDEL-PrecisionRecall.svg ${RTG_ROCPLOT_OPTIONS} --title="INDELs, model ${MODEL_TIME}"
dieIfError "Unable to generate indel Precision Recall plot."

rtg rocplot ${RTG_OUTPUT_FOLDER}/indel/non_snp_roc.tsv.gz --svg ${RTG_OUTPUT_FOLDER}/indel/INDEL-ROC.svg  ${RTG_ROCPLOT_OPTIONS} --title="INDELs, model ${MODEL_TIME}"
dieIfError "Unable to generate indel ROC plot."

# Following is currently disabled.
exit 0
# SNPs and indels, ignoring ploydy mistakes:
rtg vcfeval --baseline=${GOLD_STANDARD_VCF_GZ}  \
        -c ${VCF_OUTPUT_SORTED}.gz -o ${RTG_OUTPUT_FOLDER}/squash-ploidy --output-mode combine --template=${RTG_TEMPLATE}  \
            --evaluation-regions=${GOLD_STANDARD_CONFIDENT_REGIONS_BED_GZ} \
            --bed-regions=${BED_OBSERVED_REGIONS_OUTPUT}-sorted.bed.gz \
            --vcf-score-field=P  --sort-order=descending --squash-ploidy
dieIfError "Failed to run rtg vcfeval."

rtg rocplot ${RTG_OUTPUT_FOLDER}/squash-ploidy/non_snp_roc.tsv.gz --svg ${RTG_OUTPUT_FOLDER}/squash-ploidy/INDEL-ROC.svg
dieIfError "Unable to generate SNP ROC plot."

rtg rocplot ${RTG_OUTPUT_FOLDER}/squash-ploidy/non_snp_roc.tsv.gz -P --svg ${RTG_OUTPUT_FOLDER}/squash-ploidy/non_snp-PrecisionRecall.svg
dieIfError "Unable to generate SNP Precision Recall plot."

rtg rocplot ${RTG_OUTPUT_FOLDER}/squash-ploidy/snp_roc.tsv.gz -P --svg ${RTG_OUTPUT_FOLDER}/squash-ploidy/snp-PrecisionRecall.svg
dieIfError "Unable to generate non-SNP Precision Recall plot."
