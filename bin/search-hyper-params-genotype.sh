#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh
assertParallelInstalled

if [ "$#" -lt 2 ]; then
   echo "Argument missing. usage: number-of-runs train-genotype.sh 10g .. -n 100000 -x 100000"
   exit 1;
fi
if [ -z "${SBI_SEARCH_PARAM_CONFIG+set}" ]; then
    SBI_SEARCH_PARAM_CONFIG=config.txt

cat << EOF | cat> config.txt
--num-layers
int
3
15

--mini-batch-size
categorical
128
2048

--indel-sequence-length
int
1
20

--num-lstm-layers
int
1
5

--num-lstm-nodes-indels
int
2
20

--reduction-rate
uniform
0.3
1.3

--model-capacity
uniform
0.3
1.3

--regularization-rate
log-uniform
1E-12
1E-1

--random-seed
categorical
32434

--learning-rate
log-uniform
0.1
10

--early-stopping-num-epochs
categorical
2

--feature-mapper
categorical
org.campagnelab.dl.genotype.mappers.GenotypeMapperV35

--genomic-context-length
int
21
41

--label-smoothing-epsilon
uniform
0
0.2

--decision-threshold
uniform
0.4
0.6

--net-architecture
categorical
org.campagnelab.dl.genotype.learning.architecture.graphs.GenotypeSixDenseLayersNarrower2

EOF
    echo "SBI_SEARCH_PARAM_CONFIG not set. Using default hyper parameters. Change the variable a file with an arg-generator config file to customize the search."
fi

cat << EOF | cat>gpu.txt
0
1
2
3
EOF

echo $* >main-command.txt
NUM_GPUS=`wc -l gpu.txt|cut -d " " -f 1`

num_executions=${memory_requirement}
echo "Building caches"
arg-generator.sh 1g --config ${SBI_SEARCH_PARAM_CONFIG} --output gen-args.txt --num-commands ${num_executions}

parallel  echo `cat main-command.txt`  \
  --build-cache-then-stop \
  :::: gen-args.txt ::: \
>build-cache-commands.txt

export FORCE_PLATFORM=native
cat build-cache-commands.txt |parallel -j${NUM_GPUS} --progress --eta  --bar

echo "Training.."
unset FORCE_PLATFORM
parallel echo `cat main-command.txt` \
        --memory-cache none  --max-epochs 40 \
        :::: gen-args.txt \
>commands.txt

shuf commands.txt  |head -${num_executions} >commands-head-${num_executions}
chmod +x commands-head-${num_executions}
cat ./commands-head-${num_executions} |parallel --xapply echo :::: - ::: --gpu-device :::: gpu.txt  >all-commands.txt
cat all-commands.txt |parallel --eta --progress --bar -j${NUM_GPUS} --progress
sort -n -k 2 model-conditions.txt|tail
