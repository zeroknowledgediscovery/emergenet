#!/bin/bash

INFILE=$1

shift
count=1
while [ $# -gt 0 ]
do
    pdfjam $INFILE $1 -o page"$1".pdf
    ./geteps.sh page"$1".pdf
    mv page"$1".eps ED_Table_"$count"'.eps'
    count=$((count+1))
    shift
done   
