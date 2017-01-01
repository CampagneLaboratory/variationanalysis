#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh
assertParallelInstalled

cat << EOF | cat> learn.txt
5
EOF

cat << EOF | cat> drop.txt
0.5
1
EOF


cat << EOF | cat> reg.txt
0
1E-10
1E-7
1E-5
1E-4
1E-2
EOF

cat << EOF | cat>seed.txt
443
EOF

cat << EOF | cat>genomic-context-size.txt
21
31
41
51
61
EOF


cat << EOF | cat>is-variant-loss-weight.txt
0
5
10
50
100
1000
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

parallel echo `cat main-command.txt` --build-cache-then-stop \
  --genomic-context-length :::: genomic-context-size.txt ::: \
>build-cache-commands.txt

export FORCE_PLATFORM=native
cat build-cache-commands.txt |parallel -j${NUM_GPUS} --progress

unset FORCE_PLATFORM
parallel echo `cat main-command.txt` --regularization-rate :::: reg.txt :::  \
  --random-seed :::: seed.txt ::: --learning-rate :::: learn.txt ::: \
  --dropout-rate :::: drop.txt ::: --early-stopping-num-epochs ::: 1 ::: \
  --variant-loss-weight :::: is-variant-loss-weight.txt ::: \
  --genomic-context-length :::: genomic-context-size.txt ::: \
>commands.txt

shuf commands.txt  |head -${num_executions} >commands-head-${num_executions}
chmod +x commands-head-${num_executions}
cat ./commands-head-${num_executions} |parallel --xapply ::: echo  :::: - ::: --gpu-device :::: gpu.txt  >all-commands.txt
cat all-commands.txt |parallel -j${NUM_GPUS} --progress
sort -n -k 2 model-conditions.txt|tail
