#!/usr/bin/env bash

# PURPOSE
# generate lists of known planet host stars that TESS will observe.

sublist="knownplanet"

# sn=0
# python get_targets.py -sn $sn --get_TIC71 --get_kane_knownplanets \
#             -ts $sublist > ../results/Kane_"$sublist"_sector$sn.txt

for sn in {0..12}; do
  echo $sn ;
  python get_targets.py -sn $sn \
                        --get_TIC71 \
                        --get_kane_knownplanets \
                        -ts $sublist > ../results/Kane_"$sublist"_sector$sn.txt
done
