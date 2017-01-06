#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh
assertParallelInstalled

cat << EOF | cat> config.txt
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
1

--feature-mapper
categorical
org.campagnelab.dl.genotype.mappers.GenotypeMapperV22
org.campagnelab.dl.genotype.mappers.GenotypeMapperV23
org.campagnelab.dl.genotype.mappers.GenotypeMapperV24

--variant-loss-weight
log-uniform
0
1000

--genomic-context-length
int
21
41
EOF

cat << EOF | cat>gpu.txt
0
1
2
3
EOF




echo $* >main-command.txt
NUM_GPUS=`wc -l gpu.txt|cut -d " " -f 1`

num_executions=${memory_requirement}

. `dirname "${BASH_SOURCE[0]}"`/arg-generator.sh 1g --config config.txt --output gen-args.txt --num-commands ${num_executions}

parallel echo `cat main-command.txt` --mini-batch-size 2048 \
  --build-cache-then-stop \
  :::: gen-args.txt ::: \
>build-cache-commands.txtbuil

export FORCE_PLATFORM=native
cat build-cache-commands.txt |parallel -j${NUM_GPUS} --progress

unset FORCE_PLATFORM
parallel echo `cat main-command.txt` --mini-batch-size 2048 --memory-cache training,validation \
  :::: gen-args.txt \
>commands.txt

shuf commands.txt  |head -${num_executions} >commands-head-${num_executions}
chmod +x commands-head-${num_executions}
cat ./commands-head-${num_executions} |parallel --xapply echo :::: - ::: --gpu-device :::: gpu.txt  >all-commands.txt
cat all-commands.txt |parallel -j${NUM_GPUS} --progress
sort -n -k 2 model-conditions.txt|tail
