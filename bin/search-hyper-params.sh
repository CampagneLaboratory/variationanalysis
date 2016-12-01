#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

cat << EOF | cat> drop.txt
0.05
0.1
0.15
0.2
0.25
0.3
0.35
0.4
0.45
0.5
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
1
436
EOF

echo ${memory_requirement} $* >main-command.txt

parallel echo `cat main-command.txt` --regularization-rate :::: reg.txt :::  --random-seed :::: seed.txt ::: --dropout-rate :::: drop.txt  >commands.txt
shuf commands.txt  |head -100 >commands-head-100
chmod +x commands-head-100
./commands-head-100

