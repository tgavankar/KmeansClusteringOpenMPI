#!/bin/bash

k=(2 4 8 12)
p=(10 100 10000)

for kv in ${k[*]}
do
    for pv in ${p[*]}
    do
        printf "%s " $kv
        printf "%s\n" $pv
        python -m cProfile sequential.py -k ${kv} -u 0.0001 -i generators/input/graphk${kv}p${pv}.csv -o output/graphk${kv}p${pv}seq.csv > cProfile/graphk${kv}p${pv}seq.txt
    done
done

kv=10
pv=100000

python -m cProfile sequential.py -k ${kv} -u 0.0001 -i generators/input/graphk${kv}p${pv}.csv -o output/graphk${kv}p${pv}seq.csv > cProfile/graphk${kv}p${pv}seq.txt

for kv in ${k[*]}
do
    for pv in ${p[*]}
    do
        printf "%s " $kv
        printf "%s\n" $pv
        python -m cProfile sequential.py -t dna -k ${kv} -u 0.0001 -i generators/input/dnak${kv}p${pv}.csv -o output/dnak${kv}p${pv}seq.csv > cProfile/dnak${kv}p${pv}seq.txt
    done
done

kv=10
pv=100000

python -m cProfile sequential.py -t dna -k ${kv} -u 0.0001 -i generators/input/dnak${kv}p${pv}.csv -o output/dnak${kv}p${pv}seq.csv > cProfile/dnak${kv}p${pv}seq.txt

