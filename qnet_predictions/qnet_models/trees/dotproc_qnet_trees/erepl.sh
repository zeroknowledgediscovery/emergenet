#!/bin/bash


HELP0='From standard quasinet tree output:'
HELP1='Remove "Err" line, replace "Prob" with "frac" by default'

if [ $# -eq 0 ] ; then
    echo $HELP0
    echo $HELP1
    exit
fi

if [ $1 == '-h' ] ; then
    echo $HELP0
    echo $HELP1
    exit
fi
REPL='frac:'
LINE="index:\n"

if [ $# -gt 0 ] ; then
    REPL=$1
fi
if [ $# -gt 1 ] ; then
    REPL=$1
    LINE=$2
fi

for i in `ls *dot`
do

cat $i | grep -v "^Prob"  > proc"$i"
# remove 'x' from infront of labels
sed  -i 's/label="x\([0-9]\)/label="\1/g'  proc"$i"

j=proc"$i"

sed -i "s/Frac:/$REPL/g" $j

sed -i "s/\(label=\(\"\)\)\([0-9]\)/\1$LINE\3/g" $j

dot -Tpdf proc"$i" -o ${j/dot/pdf}
done

