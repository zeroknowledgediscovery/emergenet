#!/bin/bash

DRYRUN=0

SUBTYPE="h1n1 h1n2 h3n2 h5n1 h5n2 h5n6 h5n8 h7 h9n2"

for i in `echo $SUBTYPE`
do
    basedtaname=IXC"$i"
    PROG=" module load python; python3 animal_predictions.py $i"
    echo -n '... launching '
    T=30
    LAUNCH='/project2/ishanu/LAUNCH_UTILITY/launcher_s.sh '
    $LAUNCH -d $DRYRUN -P "$PROG"  -T $T -N 1 -C 28 -p  broadwl -J $basedtaname -M 50
done

rm *depx
