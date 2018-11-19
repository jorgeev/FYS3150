#!/bin/bash

declare -a Lattice=("40" "60" "80" "100")
t0='2.27'
ti='2.32'
step='10'
cycles='1e6'

#for L in Lattice;
#do python3 ising.py $L $t0 $ti $step $cycles; done 

for L in "${Lattice[@]}";
do python3 ising.py $L $t0 $ti $step $cycles; done

python3 latticeplot.py
