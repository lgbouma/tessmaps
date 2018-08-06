#!/usr/bin/env bash

##########
# PURPOSE
# Generate all the possible plots. You can also comment at-will if you want
# some plots, rather than others.
##########

python plot_skymaps_of_targets.py --sanity_check

python plot_skymaps_of_targets.py --all_sky

python plot_skymaps_of_targets.py --rectmap_plain -sn 0

python plot_skymaps_of_targets.py --kharchenko_2013

python plot_skymaps_of_targets.py --known_planets
