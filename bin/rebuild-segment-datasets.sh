#!/usr/bin/env bash

(cd ~/goby3; git stash; git pull; mvn install)
(cd ~/variationanalysis; git stash; git pull; ./build-all.sh -DskipTests=true)
DATE=`date +%Y-%m-%d`
PREFIX="NA12878-GIAB-chr16-chr18-chr20"
sbi-to-ssi.sh 10g -i "${PREFIX}-training.sbi" -o "${PREFIX}-training-${DATE}" -g 10 --map-features -s INDEL1 --sampling-rate 0.01 $@
sbi-to-ssi.sh 10g -i "${PREFIX}-validation.sbi" -o "${PREFIX}-validation-${DATE}" -g 10 --map-features -s INDEL1 --sampling-rate 0.01 $@
sbi-to-ssi.sh 10g -i "${PREFIX}-test.sbi" -o "${PREFIX}-validation-${DATE}" -g 10 --map-features -s INDEL1 --sampling-rate 0.01 $@

randomize-ssi.sh 10g -i "${PREFIX}-training-${DATE}" -o "${PREFIX}-random-${DATE}-train"
randomize-ssi.sh 10g -i "${PREFIX}-validation-${DATE}" -o "${PREFIX}-random-${DATE}-validation"
randomize-ssi.sh 10g -i "${PREFIX}-test-${DATE}" -o "${PREFIX}-random-${DATE}-test"