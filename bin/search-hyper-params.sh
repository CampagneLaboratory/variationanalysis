#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

cat << EOF | cat> drop.txt
0
0.3
0.4
0.5
0.6
0.7
0.8
0.9
1
EOF


cat << EOF | cat> reg.txt
0
1E-10
1E-9
1E-8
1E-7
1E-6
1E-5
1E-4
1E-3
1E-2
1E-1
EOF

cat << EOF | cat>seed.txt
2389283
443
732
EOF


echo $* >main-command.txt
cat << EOF | cat>gpu.txt
0
1
2
3
EOF

num_executions=${memory_requirement}

parallel echo `cat main-command.txt` --regularization-rate :::: reg.txt :::  --random-seed :::: seed.txt ::: --dropout-rate :::: drop.txt ::: --early-stopping-num-epochs ::: 0   >commands.txt
shuf commands.txt  |head -${num_executions} >commands-head-${num_executions}
chmod +x commands-head-${num_executions}
cat ./commands-head-${num_executions} |parallel --xapply ::: echo  :::: - ::: --gpu-device :::: gpu.txt  >all-commands.txt
cat all-commands.txt |parallel -j4
sort -n -k 2 model-conditions.txt|tail
