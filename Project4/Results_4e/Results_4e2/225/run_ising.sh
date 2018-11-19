#!/bin/bash

declare -a Lattice=("40" "60" "80" "100")
t0='2.25'
ti='2.35'
step='10'
cycles='1e7'

#for L in Lattice;
#do python3 ising.py $L $t0 $ti $step $cycles; done 

for L in "${Lattice[@]}";
do python3 ising.py $L $t0 $ti $step $cycles; done

python3 latticeplot.py
