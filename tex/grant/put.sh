#!/bin/bash

TM=`date | sed 's/ /_/g'`
sudo cp main.pdf /var/www/html/data/nrm/nrm-"$TM".pdf
