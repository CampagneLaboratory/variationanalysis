#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh


echo $* >main-command.txt
cat << EOF | cat>gpu.txt
0
1
2
3
EOF

NUM_GPUS=`wc -l gpus.txt|cut -d " " -f 1`
INPUT=${memory_requirement}
MODELS=models/*
MODEL_TIMES=`grep -v Tag model-conditions.txt | cut -d" " -f 6|awk '{print "models/"$1}' `

parallel --progress -j${NUM_GPUS} --xapply ${DLVA_HOME}/bin/predict.sh 10g -l bestAUC --records-for-auc 100000 --num-examples 100000 -f -i ${INPUT} -m ::: ${MODEL_TIMES} ::: --gpu-device :::: gpus.txt