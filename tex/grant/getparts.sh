#!/bin/bash

FILE=$1

pdfjam $FILE 1-5 --letterpaper -o narrative.pdf
pdfjam $FILE 6 --letterpaper -o abbreviations.pdf
pdfjam $FILE 7-8 --letterpaper -o datamanagement.pdf
pdfjam $FILE 9 --letterpaper -o technicalabstract.pdf
pdfjam $FILE 10 --letterpaper -o layabstract.pdf
pdfjam $FILE 11-13 --letterpaper -o sow.pdf
pdfjam $FILE 14 --letterpaper -o impact.pdf
pdfjam $FILE 15 --letterpaper -o militaryhealth.pdf
pdfjam $FILE 16-20 --letterpaper -o facilities.pdf
pdfjam $FILE 21 --letterpaper -o patentspubs.pdf
pdfjam $FILE 24-25 --letterpaper -o IP.pdf
pdfjam $FILE 26 --letterpaper -o datasharing.pdf
pdfjam $FILE 27- --letterpaper -o refrences.pdf




