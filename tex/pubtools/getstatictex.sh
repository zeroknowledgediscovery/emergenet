#!/bin/bash

FILE=$1
cp $FILE ${FILE/.tex/static.tex}

grep newlabel ${FILE/tex/aux} | sed 's/\\newlabel//g' | awk -F"\}\{" '{print $1,$2}' | sed 's/{//g' | grep -v \# > _tmp

awk 'NF==2{print}' _tmp > _tmp1
L=`wc -l _tmp1 | awk '{print $1}'`

l=1

while [ $l -le $L ]
do
    ln=`sed -n "$l"p _tmp1`
    a=`echo $ln | awk '{print $1}'`
    b=`echo $ln | awk '{print $2}'`

    sed -i "s/$a/$b/g" ${FILE/.tex/static.tex}
    
    l=$((l+1))
done

sed -i 's/\\ref{\([0-9]*\)}/\1/g'  ${FILE/.tex/static.tex}
rm _tmp*
