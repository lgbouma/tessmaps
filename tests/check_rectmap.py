'''
run this script from anywhere on your system. change `savdir` below!
'''
import numpy as np
from tessmaps import tessmaps as tm

from astropy.coordinates import SkyCoord
coords = SkyCoord([42, 42.42, 241.04, 225], [-68, -60, -22.58, -28], unit='deg')

names = np.array(['hi','friend','obj','obj'])
is_transiting = np.array([True, True, True, True])
savname = 'temp.png'
savdir = '/Users/luke/' # fix this to your preferred directory!
title = 'my map'
sector_number = 10

tm.make_rect_map(sector_number, coords, names=names,
                 annotate_bools=is_transiting, title=title,
                 bkgnd_cmap='Blues', savname=savname, savdir=savdir)

print('if you got here w/out errors, it worked!')
