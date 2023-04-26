#!/bin/bash

FILE=$1
mkdir BOX >& /dev/null

pdfjam $FILE 1-5 --letterpaper -o X_narrative.pdf
pdfjam $FILE 6 --letterpaper -o X_abbreviations.pdf
pdfjam $FILE 7-8 --letterpaper -o X_datamanagement.pdf
pdfjam $FILE 9 --letterpaper -o X_technicalabstract.pdf
pdfjam $FILE 10 --letterpaper -o X_layabstract.pdf
pdfjam $FILE 11-13 --letterpaper -o X_sow.pdf
pdfjam $FILE 14 --letterpaper -o X_impact.pdf
pdfjam $FILE 15 --letterpaper -o X_militaryhealth.pdf
pdfjam $FILE 16-20 --letterpaper -o X_facilities.pdf
pdfjam $FILE 21 --letterpaper -o X_patentspubs.pdf
pdfjam $FILE 22-23 --letterpaper -o X_LOS1.pdf
pdfjam $FILE 24 --letterpaper -o X_LOS2.pdf
pdfjam $FILE 25-26 --letterpaper -o X_IP.pdf
pdfjam $FILE 27-28 --letterpaper -o X_datasharing.pdf
pdfjam $FILE 29- --letterpaper -o X_references.pdf


mv X_*pdf BOX
