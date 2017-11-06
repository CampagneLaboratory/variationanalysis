#!/usr/bin/env bash

(cd ~/goby3; git stash; git pull; mvn install)
(cd ~/variationanalysis; git stash; git pull; ./build-all.sh -DskipTests=true)
DATE=`date +%Y-%m-%d`
sbi-to-ssi.sh 10g -i out-part-11.sbi -o training-${DATE} -g 10 --map-features -s INDEL1 --sampling-rate 0.01
sbi-to-ssi.sh 10g -i out-part-41.sbi -o validation-${DATE} -g 10 --map-features -s INDEL1 --sampling-rate 0.01
