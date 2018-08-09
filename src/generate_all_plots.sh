#!/usr/bin/env bash

##########
# PURPOSE
# Generate all the possible plots. Comment at-will if you want some plots,
# rather than others.
##########

python plot_skymaps_of_targets.py --sanity_check

python plot_skymaps_of_targets.py --all_sky

for sn in {0..12}; do
  python plot_skymaps_of_targets.py --rectmap_plain -sn $sn
done

python plot_skymaps_of_targets.py --known_planets
 
python plot_skymaps_of_targets.py --kharchenko_2013

