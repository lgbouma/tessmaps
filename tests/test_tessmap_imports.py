'''
if things break here, it's because the package setup did not work
'''

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

from tessmaps import tessmaps as tm

# check: did we successfully export the TIC7.1 cut data files?
elon, elat, totsec = tm._get_TIC_coords_count(5)

# and did we successfully export the other important source files?
coords = SkyCoord(lon=np.array([0,0])*u.deg, lat=np.array([-90,-90])*u.deg,
                  frame='barycentrictrueecliptic')
coords = coords.icrs
tm.make_rect_map(4, coords,title='hi',savname='temp.png',savdir='')

from tessmaps import get_time_on_silicon as gts
df = gts.get_time_on_silicon(coords)

from tessmaps import get_targets as gt

print('if you got here w/out errors, it worked!')
