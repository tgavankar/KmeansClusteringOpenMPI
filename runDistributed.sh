#!/bin/bash

np=(2 4 8 12)
k=(2 4 8 12)
p=(10 100 10000)

for npv in ${np[*]}
do
    for kv in ${k[*]}
    do
        for pv in ${p[*]}
        do
            printf "%s " $kv
            printf "%s\n" $pv
            mpirun --machinefile machines -np ${npv} python distributed.py -k ${kv} -u 0.0001 -i generators/input/graphk${kv}p${pv}.csv -o output/graphk${kv}p${pv}mpi.csv > cProfile/graphk${kv}p${pv}np${npv}mpi.txt
        done
    done


kv=10
pv=100000

mpirun --machinefile machines -np ${npv} python distributed.py -k ${kv} -u 0.0001 -i generators/input/graphk${kv}p${pv}.csv -o output/graphk${kv}p${pv}mpi.csv > cProfile/graphk${kv}p${pv}np${npv}mpi.txt


for kv in ${k[*]}
do
    for pv in ${p[*]}
    do
        printf "%s " $kv
        printf "%s\n" $pv
        mpirun --machinefile machines -np ${npv} python distributed.py -t dna -k ${kv} -u 0.0001 -i generators/input/dnak${kv}p${pv}.csv -o output/dnak${kv}p${pv}mpi.csv > cProfile/dnak${kv}p${pv}np${npv}mpi.txt
    done
done

kv=10
pv=100000

mpirun --machinefile machines -np ${npv} python distributed.py -t dna -k ${kv} -u 0.0001 -i generators/input/dnak${kv}p${pv}.csv -o output/dnak${kv}p${pv}mpi.csv > cProfile/dnak${kv}p${pv}np${npv}mpi.txt

done
