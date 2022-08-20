#!/bin/bash

echo -n "nature med word count: intro+result+discussion: "

texcount main_new.tex | head -n 20 | tail -n 8 | awk -F+ 'BEGIN{s=0}{s=s+$1}END{print s}'
