#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh
assertParallelInstalled

if [ "$#" -lt 1 ]; then
   echo "Argument missing. usage: number-of-runs train-autoencoder.sh  "
   echo "For instance:  ~/variationanalysis/bin/search-hyper-params-genotypetensors.sh 10 --problem genotyping:CNG-NA12878-realigned-2018-01-30 --mode supervised_genotypes --num-workers 0 --autoencoder-type 1 "
   echo "Monitor performance with tail --lines 1 best-perfs-*|sort -k 9 -n"
   exit 1;
fi
if [ -z "${SBI_SEARCH_PARAM_CONFIG+set}" ]; then
    SBI_SEARCH_PARAM_CONFIG=search-config.txt

cat << EOF | cat> search-config.txt
--lr
log-uniform
1e-01
1e-02

--num-layers
int
1
5

--L2
log-uniform
1E-1
1E-20

--dropout-probability
uniform
0
0.2

--mode
categorical
supervised_genotypes

--mini-batch-size
categorical
128

--encoded-size
int
4
128

--indel-weight-factor
log-uniform
1
100

--max-epochs
categorical
20

-n
categorical
100000

-x
categorical
100000

EOF
    echo "SBI_SEARCH_PARAM_CONFIG not set. Using default hyper parameters. Change the variable a file with an arg-generator config file to customize the search."
fi

cat << EOF | cat>gpu.txt
0
1
2
3
4
5
6
EOF

echo $* >main-command.txt
NUM_GPUS=`wc -l gpu.txt|cut -d " " -f 1`

num_executions=${memory_requirement}

arg-generator.sh 1g --config ${SBI_SEARCH_PARAM_CONFIG} --output gen-args.txt --num-commands ${num_executions}

echo "Training.."

parallel echo `cat main-command.txt` \
        --max-epochs 1000 \
        :::: gen-args.txt \
>commands.txt

shuf commands.txt  |head -${num_executions} >commands-head-${num_executions}
chmod +x commands-head-${num_executions}
cat ./commands-head-${num_executions}  |parallel --trim lr --xapply echo  run-on-gpu.sh :::: gpu.txt ::::  -   >all-commands.txt
cat all-commands.txt |parallel --line-buffer --eta --progress --bar -j${NUM_GPUS} --progress
sort -n -k 10 best-perfs-*.tsv |tail -10
COMMANDS=./commands-head-${num_executions}-${RANDOM}.txt
cp all-commands.txt ${COMMANDS}
echo "Commands have been saved in ${COMMANDS}"


