#!/bin/bash

NAME='main.pdf'

if [ $# -gt 0 ] ; then
    NAME=$1
fi
cp $NAME IPF.pdf

sftp  ishanu@webspace.uchicago.edu:zed/data/pub_drafts_/PAI/ <<< $'put -r IPF.pdf'
