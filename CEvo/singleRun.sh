#!/bin/bash
for j in $(seq $1 $1); do
let "x=$j"
let "i=$x * $x"

export x

cat submit.sh.template | sbatch --ntasks=$i
done