#!/bin/bash

FILE=main.tex

if [ $# -gt 0 ]; then
   FILE=$1
fi

sed 's/\(ref{\)/\n\1/g' $FILE | sed -n 's/.*ref{\([^}]*\).*/\1/gp' | awk '!count[$0]++'
