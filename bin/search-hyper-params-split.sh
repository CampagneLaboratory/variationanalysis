#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh
assertParallelInstalled

if [ "$#" -lt 1 ]; then
   echo "Argument missing. usage: number-of-runs train-with-ureg3-*.sh  "
   exit 1;
fi
if [ -z "${SBI_SEARCH_PARAM_CONFIG+set}" ]; then
    SBI_SEARCH_PARAM_CONFIG=config.txt

cat << EOF | cat> config.txt
--lr
log-uniform
1e-01
1e-05

--lr-patience
int
3
10

--factor
log-uniform
10
0.01

--L2
log-uniform
1E-2
1E-12

--mode
categorical
split

--constant-learning-rates
flag

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
        --max-epochs 200 --abort-when-failed-to-improve 20 \
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


