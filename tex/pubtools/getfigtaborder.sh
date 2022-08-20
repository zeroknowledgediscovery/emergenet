#!/bin/bash

FILE=main.aux

if [ $# -gt 0 ] ; then
    FILE=$1
fi
grep newlabel $FILE | sed -n  's/}{/ /gp' | awk '{print $1,$2,$3}' | sed 's/\\newlabel{//g' | grep -v \#  | sed 's/{//g'
