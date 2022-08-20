#!/bin/bash

if [ $# -lt 2 ] ; then
SI='SI.pdf'
else
SI=$2
fi

pdfjam --paper letterpaper --rotateoversize false $1 $SI -o authorpdf.pdf
