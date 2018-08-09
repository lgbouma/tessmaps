'''
Some example uses:

from tessmaps import tessmaps as tm

df = tm.get_time_on_silicon(coords)

tm.make_rect_map(sector_number, coords, names=names,
                 annotate_bools=is_transiting, title=title,
                 bkgnd_cmap='Blues', savname=savname)

tm.make_sector_list(sector_number, coords, names=names,
                    annotate_bools=is_transiting, title=title,
                    bkgnd_cmap='Blues', savname=savname)
'''

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

from . import tessmaps, get_time_on_silicon, plot_skymaps_of_targets
