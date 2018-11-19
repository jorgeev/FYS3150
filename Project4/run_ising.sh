#!/bin/bash

declare -a Lattice=("40" "60" "80" "100")
t0='2.0'
ti='2.4'
step='16'
cycles='1e7'

#for L in Lattice;
#do python3 ising.py $L $t0 $ti $step $cycles; done 

for L in "${Lattice[@]}";
do python3 ising.py $L $t0 $ti $step $cycles; done
