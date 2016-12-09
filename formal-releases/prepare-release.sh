#!/usr/bin/env bash

DL_VERSION=`cat ../VERSION.txt`
WORKDIR=`pwd`
rm -rf release-dlvariation_$DL_VERSION
mkdir -p release-dlvariation_$DL_VERSION
cd release-dlvariation_$DL_VERSION
cp ../../VERSION.txt .
cp -r ../../bin .
mkdir -p gpus/target
cp ../../gpus/target/gpus-$DL_VERSION.jar gpus/target/
mkdir -p somatic/target
cp ../../somatic/target/somatic-$DL_VERSION-bin.jar somatic/target/
cd ${WORKDIR}