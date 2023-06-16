#!/bin/bash
TAG=$*

for tag in `echo $TAG`
	   do
	       NAME='enetshap'"$tag"
	       sbatch --job-name=$NAME --output=OUT$NAME --error=ERR$NAME --time=20:00:00  --mem=100G --cpus-per-task=28 --ntasks=1 --wrap="module load gcc/11.3.0; module load python/3.10.5;python3.10 ./shapanalysis_single_sequence.py --sequence_name $tag"
done

#./runsbatch.sh north_h1n1_ha north_h3n2_ha south_h1n1_ha south_h3n2_ha
