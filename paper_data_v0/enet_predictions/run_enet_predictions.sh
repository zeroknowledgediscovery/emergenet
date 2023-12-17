#!/bin/bash

DRYRUN=0

HEMISPHERE="north south"
SUBTYPE="h1n1 h3n2"

for i in `echo $HEMISPHERE`
do
    for j in `echo $SUBTYPE`
    do
        for ((k=0;k<=20;k++)); 
        do
            basedtaname=IXC"$i$j$k"
            PROG=" module load python; python3 enet_predictions.py $i $j $k"
            echo -n '... launching '
            T=30
            LAUNCH='/project2/ishanu/LAUNCH_UTILITY/launcher_s.sh '
            $LAUNCH -d $DRYRUN -P "$PROG"  -T $T -N 1 -C 28 -p  broadwl -J $basedtaname -M 50
        done
    done
done

rm *depx
