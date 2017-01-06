#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/common.sh

if [ "$#" -ne 3 ]; then
   echo "Usage: evaluate-genotypes.sh model-directory model-prefix test-set.sbi."
   echo "The env variables GOLD_STANDARD_VCF_GZ and GOLD_STANDARD_CONFIDENT_REGIONS_BED_GZ can be used to change the VCF and confident region bed."
   echo "The first run downloads these files from the Genome in a Bottle for sample NA12878 when the variables are not defined."
   exit 1;
fi
MODEL_DIR=$1
MODEL_PREFIX=$2
DATASET_SBI=$3

if [ -e configure.sh ]; then
 echo "Loading configure.sh"
 source configure.sh
fi
if [ ! -e hg19.sdf ]; then
 echo "You must install hg19.sdf in the current directory. See rtg downloads at http://www.realtimegenomics.com/news/pre-formatted-reference-datasets/"
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

if [ -z "${GOLD_STANDARD_VCF_GZ+set}" ]; then
    echo "Downloading Gold standard Genome in a Bottle VCF"
    wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.1/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf_phased.vcf.gz
    mv HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf_phased.vcf.gz GIAB-NA12878-confident.vcf.gz
    # add "chr prefix:"
    gzip -c -d GIAB-NA12878-confident.vcf.gz|awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' >GIAB-NA12878-confident-chr.vcf
    # remove non-SNPs:
    cat GIAB-NA12878-confident-chr.vcf|awk '{if($0 !~ /^#/) { if (length($4)==1 && length($5)==1) print $0;}  else {print $0}}' >GIAB-NA12878-confident-chr-snps.vcf
    bgzip GIAB-NA12878-confident-chr-snps.vcf
    tabix GIAB-NA12878-confident-chr-snps.vcf.gz
    rm GIAB-NA12878-confident-chr.vcf

    GOLD_STANDARD_VCF_GZ="GIAB-NA12878-confident-chr-snps.vcf.gz"
    echo 'export GOLD_STANDARD_VCF_GZ="GIAB-NA12878-confident-chr-snps.vcf.gz"' >>configure.sh
    echo "Gold standard VCF downloaded for NA12878 (SNPs) and named in configure.sh. Edit GOLD_STANDARD_VCF_GZ to switch to a different gold-standard validation VCF."
fi


if [ -z "${GOLD_STANDARD_CONFIDENT_REGIONS_BED_GZ+set}" ]; then
    echo "Downloading Gold standard Genome in a Bottle Confident Regions (bed)"
    wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.1/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf.bed
    mv HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.1_highconf.bed GIAB-NA12878-confident-regions.bed
    # add "chr prefix:"
    cat  GIAB-NA12878-confident-regions.bed |awk '{print "chr"$1"\t"$2"\t"$3}' >GIAB-NA12878-confident-regions-chr.bed
    bgzip GIAB-NA12878-confident-regions-chr.bed
    tabix GIAB-NA12878-confident-regions-chr.bed.gz
    rm GIAB-NA12878-confident-regions.bed

    GOLD_STANDARD_CONFIDENT_REGIONS_BED_GZ="GIAB-NA12878-confident-regions-chr.bed.gz"
    echo 'export GOLD_STANDARD_CONFIDENT_REGIONS_BED_GZ="GIAB-NA12878-confident-regions-chr.bed.gz"' >>configure.sh
    echo "Gold standard confident regions downloaded for NA12878  and named in configure.sh. Edit GOLD_STANDARD_CONFIDENT_REGIONS_BED_GZ to switch to a different gold-standard confident region bed file."
fi

if [ ! -e "${DATASET_SBI}" ]; then
        echo "The test set was not found: ${DATASET_SBI}  "
        exit 1;
fi

MODEL_TIME=`basename ${MODEL_DIR}`

echo "Running predict-genotypes to create VCF and observed region bed.."
predict-genotypes.sh 10g -m ${MODEL_DIR} -l ${MODEL_PREFIX} -f -i ${DATASET_SBI} \
    --format VCF -n 10000
dieIfError "Failed to predict dataset with model ${MODEL_DIR}/."
echo "Evaluation with rtg vcfeval starting.."

VCF_OUTPUT=`ls -1tr ${MODEL_TIME}-${MODEL_PREFIX}-*.vcf|tail -1`
BED_OBSERVED_REGIONS_OUTPUT=`ls -1tr ${MODEL_TIME}-${MODEL_PREFIX}-*-observed-regions.bed |tail -1`

VCF_OUTPUT_SORTED=`basename ${VCF_OUTPUT} .vcf`-sorted.vcf
cat ${VCF_OUTPUT} | vcf-sort > ${VCF_OUTPUT_SORTED}
dieIfError "Unable to sort prediction VCF."

bgzip -f ${VCF_OUTPUT_SORTED}
tabix ${VCF_OUTPUT_SORTED}.gz

bedtools sort -i ${BED_OBSERVED_REGIONS_OUTPUT} > ${BED_OBSERVED_REGIONS_OUTPUT}-sorted.bed
bgzip -f ${BED_OBSERVED_REGIONS_OUTPUT}-sorted.bed
tabix ${BED_OBSERVED_REGIONS_OUTPUT}-sorted.bed.gz

RTG_OUTPUT_FOLDER=output-${RANDOM}

rtg vcfeval --baseline=${GOLD_STANDARD_VCF_GZ}  \
        -c ${VCF_OUTPUT_SORTED}.gz -o ${RTG_OUTPUT_FOLDER} --template=hg19.sdf \
            --evaluation-regions=${GOLD_STANDARD_CONFIDENT_REGIONS_BED_GZ} \
            --bed-regions=${BED_OBSERVED_REGIONS_OUTPUT}-sorted.bed.gz \
            --vcf-score-field=P  --sort-order=descending
dieIfError "Failed to run rtg vcfeval."
echo "See rtg vcfeval detailed output in ${RTG_OUTPUT_FOLDER}"

rtg rocplot ${RTG_OUTPUT_FOLDER}/snp_roc.tsv.gz --svg ${RTG_OUTPUT_FOLDER}/SNP-ROC.svg
                                         dieIfError "Unable to generate ROC plot."

rtg rocplot ${RTG_OUTPUT_FOLDER}/snp_roc.tsv.gz -P --svg ${RTG_OUTPUT_FOLDER}/SNP-PrecisionRecall.svg
dieIfError "Unable to generate Precision Recall plot."