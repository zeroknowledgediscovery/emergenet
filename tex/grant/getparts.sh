#!/bin/bash

FILE=$1

pdfjam $FILE 1 --letterpaper -o summary.pdf
pdfjam $FILE 2- --letterpaper -o narrative.pdf
