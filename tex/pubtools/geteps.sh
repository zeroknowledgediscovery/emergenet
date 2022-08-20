#!/bin/bash
input=$1
gs -q -dNOCACHE -dNOPAUSE -dBATCH -dSAFER -sDEVICE=eps2write -sOutputFile=${input/pdf/eps} $input
