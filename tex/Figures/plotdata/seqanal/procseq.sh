#!/bin/bash
RECFILE=../h1n1humanHA_Northern_WHO_with_extra_months.csv
RESPREFIX='res.fasta'

if [ $# -gt 0 ] ; then
    RECFILE=$1
fi
if [ $# -gt 1 ] ; then
    RESPREFIX=$2
fi



awk -F, '{print $1,$3,$5,$8}' $RECFILE > tmp

LINES=`wc -l tmp | awk '{print $1}'`
line=2

while [ $line -le $LINES ]
do
    LINE=`sed -n "$line"p tmp`
    YR=`echo $LINE | cut -d " " -f 1 | sed 's/\_/-/g'`
    DM=`echo $LINE | cut -d " " -f 3`
    WH=`echo $LINE | cut -d " " -f 2`
    QN=`echo $LINE | cut -d " " -f 4`

    echo $YR

    RESFILE="$YR""$RESPREFIX"

    echo '>WHO:' > $RESFILE
    echo $WH >> $RESFILE
    echo '' >> $RESFILE

    echo '>DOM:' >> $RESFILE
    echo $DM >> $RESFILE
    echo '' >> $RESFILE

    echo '>QNT:' >> $RESFILE
    echo $QN >> $RESFILE
    echo '' >> $RESFILE

    
    line=$((line+1))
done
