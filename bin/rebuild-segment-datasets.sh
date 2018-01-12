#!/usr/bin/env bash
set -x

# (cd ~/goby3; git stash; git pull; mvn install)
# (cd ~/variationanalysis; git stash; git pull; ./build-all.sh -DskipTests=true)
DATE=`date +%Y-%m-%d`
source ./configure.sh

PREFIX="NA12878-GIAB"
GENOMIC_CONTEXT_LENGTH="--genomic-context-length 1"
VARMAP="/scratchLocal/joc2080/reference_varmaps/NA12878-GIAB-gold.varmap"

simulate-sbi.sh 10g  -i ${VARMAP} \
    --include chr4 --include chr6 --include chr8  -o NA12878-GIAB-training.sbi \
    --genome ${SBI_GENOME} ${GENOMIC_CONTEXT_LENGTH}

simulate-sbi.sh 10g  -i /scratchLocal/joc2080/reference_varmaps/NA12878-GIAB-gold.varmap \
    --include chr16  -o NA12878-GIAB-validation.sbi  \
    --genome ${SBI_GENOME} ${GENOMIC_CONTEXT_LENGTH}

simulate-sbi.sh 10g  -i /scratchLocal/joc2080/reference_varmaps/NA12878-GIAB-gold.varmap \
    --include chr20  -o NA12878-GIAB-test.sbi  \
    --genome ${SBI_GENOME} ${GENOMIC_CONTEXT_LENGTH}

OPTIONS="-g 10 --map-features -s INDEL1 --sampling-rate 0.01 ${GENOMIC_CONTEXT_LENGTH} "
sbi-to-ssi.sh 10g -i "${PREFIX}-training.sbi" -o "${PREFIX}-training-${DATE}"  ${OPTIONS} $@
sbi-to-ssi.sh 10g -i "${PREFIX}-validation.sbi" -o "${PREFIX}-validation-${DATE}" ${OPTIONS} $@
sbi-to-ssi.sh 10g -i "${PREFIX}-test.sbi" -o "${PREFIX}-validation-${DATE}"  ${OPTIONS} $@

randomize-ssi.sh 10g -i "${PREFIX}-training-${DATE}" -o "${PREFIX}-random-${DATE}-train"
randomize-ssi.sh 10g -i "${PREFIX}-validation-${DATE}" -o "${PREFIX}-random-${DATE}-validation"
randomize-ssi.sh 10g -i "${PREFIX}-test-${DATE}" -o "${PREFIX}-random-${DATE}-test"