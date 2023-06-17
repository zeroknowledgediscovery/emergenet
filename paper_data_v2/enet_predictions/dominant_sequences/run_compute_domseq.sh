#!/bin/bash

DRYRUN=0

HEMISPHERE="north south"
SUBTYPE="h1n1 h3n2"
SEGMENT="ha na"

for i in `echo $HEMISPHERE`
do
    for j in `echo $SUBTYPE`
    do
        for k in `echo $SEGMENT`
        do
            for ((l=0;l<=19;l++)); 
            do
                basedtaname=IXC"$i$j$k$l"
                PROG=" module load python/anaconda-2021.05; python3 compute_domseq.py $i $j $k $l"
                echo -n '... launching '
                T=30
                LAUNCH='/project2/ishanu/LAUNCH_UTILITY/launcher_s.sh '
                $LAUNCH -d $DRYRUN -P "$PROG"  -T $T -N 1 -C 28 -p  broadwl -J $basedtaname -M 50
            done
        done
    done
done

rm *depx
