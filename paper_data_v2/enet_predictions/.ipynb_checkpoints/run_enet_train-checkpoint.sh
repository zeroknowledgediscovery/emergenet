#!/bin/bash

DRYRUN=0

HEMISPHERE="north south"
SUBTYPE='h1n1 h3n2'
YEAR="02_03 03_04 04_05 05_06 06_07 07_08 08_09 09_10 10_11 11_12 12_13 13_14 14_15 15_16 16_17 17_18 18_19 19_20 20_21 21_22 22_23"

for i in `echo $HEMISPHERE`
do
    for j in `echo $SUBTYPE`
    do
        for k in `echo $YEAR`
        do
            basedtaname=IXC"$i$j$k"
            PROG=" module load python; python3 enet_train.py $i $j $k"
            echo -n '... launching '
            T=30
            LAUNCH='/project2/ishanu/LAUNCH_UTILITY/launcher_s.sh '
            $LAUNCH -d $DRYRUN -P "$PROG"  -T $T -N 1 -C 28 -p  broadwl -J $basedtaname -M 50
        done
    done
done

rm *depx