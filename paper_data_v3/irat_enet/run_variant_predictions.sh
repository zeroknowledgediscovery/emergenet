#!/bin/bash

DRYRUN=0

# Number of months prior to emergence to compute risk
MONTHS="0 1 3 6 12"

for ((i=0;i<=18;i++)); 
do
    for j in `echo $MONTHS`
    do
        basedtaname=IXC"$i$j"
        PROG=" module load python/anaconda-2021.05; python3 variant_predictions.py $i $j"
        echo -n '... launching '
        T=30
        LAUNCH='/project2/ishanu/LAUNCH_UTILITY/launcher_s.sh '
        $LAUNCH -d $DRYRUN -P "$PROG"  -T $T -N 1 -C 28 -p  broadwl -J $basedtaname -M 50
    done
done

rm *depx
