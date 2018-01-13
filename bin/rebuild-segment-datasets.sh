#!/usr/bin/env bash
set -x

# (cd ~/goby3; git stash; git pull; mvn install)
# (cd ~/variationanalysis; git stash; git pull; ./build-all.sh -DskipTests=true)
DATE=`date +%Y-%m-%d`
if [ -e configure.sh ]; then
 echo "Loading configure.sh"
 source configure.sh
fi

if [ -z "${PREFIX+set}" ]; then
    PREFIX="NA12878-GIAB"
    echo "PREFIX set to ${PREFIX}. Change the variable to switch the naming of the produced dataset."
fi

if [ -z "${GENOMIC_CONTEXT_LENGTH+set}" ]; then
    GENOMIC_CONTEXT_LENGTH="--genomic-context-length 1"
    echo "GENOMIC_CONTEXT_LENGTH set to ${GENOMIC_CONTEXT_LENGTH}. Change the variable to switch the context length."
fi

if [ -z "${VARMAP+set}" ]; then
    VARMAP="/scratchLocal/joc2080/reference_varmaps/NA12878-GIAB-gold.varmap"
    echo "VARMAP set to ${VARMAP}. Change the variable to switch the location of the varmap with true calls."
fi


simulate-sbi.sh 10g  -i ${VARMAP} \
    --include chr4 --include chr6 --include chr8  \
    -o "${PREFIX}-training.sbi" \
    --genome ${SBI_GENOME} ${GENOMIC_CONTEXT_LENGTH}

simulate-sbi.sh 10g  -i ${VARMAP}  \
    --include chr16  -o "${PREFIX}-validation.sbi"  \
    --genome ${SBI_GENOME} ${GENOMIC_CONTEXT_LENGTH}

simulate-sbi.sh 10g  -i ${VARMAP}  \
    --include chr20  -o "${PREFIX}-test.sbi"  \
    --genome ${SBI_GENOME} ${GENOMIC_CONTEXT_LENGTH}

OPTIONS="-g 10 --map-features -s INDEL1 --sampling-rate 0.01 ${GENOMIC_CONTEXT_LENGTH} "
sbi-to-ssi.sh 10g -i "${PREFIX}-training.sbi" -o "${PREFIX}-training-${DATE}"  ${OPTIONS} $@
sbi-to-ssi.sh 10g -i "${PREFIX}-validation.sbi" -o "${PREFIX}-validation-${DATE}" ${OPTIONS} $@
sbi-to-ssi.sh 10g -i "${PREFIX}-test.sbi" -o "${PREFIX}-test-${DATE}"  ${OPTIONS} $@

randomize-ssi.sh 10g -i "${PREFIX}-training-${DATE}" -o "${PREFIX}-random-${DATE}-train"
randomize-ssi.sh 10g -i "${PREFIX}-validation-${DATE}" -o "${PREFIX}-random-${DATE}-validation"
randomize-ssi.sh 10g -i "${PREFIX}-test-${DATE}" -o "${PREFIX}-random-${DATE}-test"

sbi-to-ssi.sh 10g -i "${PREFIX}-validation.sbi" -o "${PREFIX}-validation-2018-01-12" -g 10 --map-features -s INDEL1 --sampling-rate 0.01 --genomic-context-length 21
sbi-to-ssi.sh 10g -i "${PREFIX}-test.sbi" -o "${PREFIX}-test-2018-01-12" -g 10 --map-features -s INDEL1 --sampling-rate 0.01 --genomic-context-length 21

randomize-ssi.sh 10g -i "${PREFIX}-validation-2018-01-12" -o "${PREFIX}-random-2018-01-12-validation"

randomize-ssi.sh 10g -i "${PREFIX}-test-2018-01-12" -o "${PREFIX}-random-2018-01-12-test"
