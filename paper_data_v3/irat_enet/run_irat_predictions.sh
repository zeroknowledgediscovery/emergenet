#!/bin/bash

DRYRUN=0


for ((i=0;i<=23;i++)); 
do
    basedtaname=IXC"$i"
    PROG=" module load python/anaconda-2021.05; python3 irat_predictions.py $i"
    echo -n '... launching '
    T=30
    LAUNCH='/project2/ishanu/LAUNCH_UTILITY/launcher_s.sh '
    $LAUNCH -d $DRYRUN -P "$PROG"  -T $T -N 1 -C 28 -p  broadwl -J $basedtaname -M 50
done

rm *depx
