#!/usr/bin/env bash

DL_VERSION=`cat ../VERSION.txt`
WORKDIR=`pwd`
rm -rf release-dlvariation_$DL_VERSION*
mkdir -p release-dlvariation_$DL_VERSION
cd release-dlvariation_$DL_VERSION
cp ../../VERSION.txt .
cp ../../PROFILE.txt .
cp ../../LICENSE.md .
cp ../../README.md .
cp ../../SOMATIC-TUTORIAL.md . 
cp -r ../../bin .
mkdir -p gpus/target
cp ../../gpus/target/gpus-$DL_VERSION.jar gpus/target/
mkdir -p somatic/target
cp ../../somatic/target/somatic-$DL_VERSION-bin.jar somatic/target/
cd ${WORKDIR}
zip -r release-dlvariation_1.1.1.zip release-dlvariation_1.1.1/
