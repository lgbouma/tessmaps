#!/usr/bin/env bash

# PURPOSE
# generate lists of known planet host stars that TESS will observe.

sublist="knownplanet"

for sn in {0..12}; do
  python get_targets.py -sn $sn --get_TIC71 -ts $sublist > ../results/"$sublist"_sector$sn.txt
done
