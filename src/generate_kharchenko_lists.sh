#!/usr/bin/env bash

# PURPOSE
# generate text lists of the best kharchenko+ 2013 clusters for each sector in
# the south.

for sn in {0..12}; do
  echo $sn
	python get_targets.py --sector_number $sn \
                        --get_kharchenko_2013 \
                        > ../results/kharchenko_sector$sn.txt
done
