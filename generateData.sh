cd generators

k=(2 4 8 12 25)
p=(100 1000 5000 10000)

for kv in ${k[*]}
do
    for pv in ${p[*]}
    do
        printf "%s " $kv
        printf "%s\n" $pv
        python graphPoints.py -c ${kv} -p ${pv} -o input/graphk${kv}p${pv}.csv
    done
done

for kv in ${k[*]}
do
    for pv in ${p[*]}
    do
        printf "%s " $kv
        printf "%s\n" $pv
        python dnaStrands.py -c ${kv} -p ${pv} -o input/dnak${kv}p${pv}.csv
    done
done







