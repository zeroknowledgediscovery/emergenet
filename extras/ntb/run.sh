#!/bin/bash
TAG=$1
python3.10 ./shapanalysis.py --hemi $TAG --SAMPLE_FRACTION .2
