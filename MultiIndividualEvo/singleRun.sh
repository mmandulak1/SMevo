#!/bin/bash
for j in $(seq $1 $1); do
let "x=$2"
let "i=$j"
let "k=$i"
export x
export k

cat submit.sh.template | sbatch --ntasks=$i
done
