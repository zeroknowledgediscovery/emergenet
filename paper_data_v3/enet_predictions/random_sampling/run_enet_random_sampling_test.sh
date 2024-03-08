#!/bin/bash

DRYRUN=0
SUBTYPE="h1n1 h3n2"
YEAR="15_16 16_17 17_18 18_19 19_20"

for i in `echo $SUBTYPE`
do
    for j in `echo $YEAR`
    do
        basedtaname=IXC"$i$j"
        PROG=" module load python/anaconda-2021.05; python3 enet_random_sampling_test.py $i $j"
        echo -n '... launching '
        T=30
        LAUNCH='/project2/ishanu/LAUNCH_UTILITY/launcher_s.sh '
        $LAUNCH -d $DRYRUN -P "$PROG"  -T $T -N 1 -C 28 -p  broadwl -J $basedtaname -M 50
    done
done

rm *depx
