Overview
----------

This repo contains some tools for understanding what TESS is looking at.

Useful pull requests are gratefully accepted. The focal plane geometry model is
not 100% accurate, but it's pretty good. 

USAGE
----------
To see if some coordinates are observed:
```
from astropy.coordinates import SkyCoord
from tessmaps import get_time_on_silicon as gts
coords = SkyCoord(['124.532, -68.313'], ['42.42, -42.42'], unit='deg')
df = gts.get_time_on_silicon(coords)
```
where `df` is a pandas DataFrame that tells you in which sector the observation
occurs, using the pointings given by
[MIT's TESS Science Office](https://tess.mit.edu/tess-pointing) and the focal
plane geometry used by
[Sullivan et al (2015)](https://arxiv.org/abs/1506.03845) and
[Bouma et al (2017)](https://arxiv.org/abs/1705.08891).

A separate module accepts sector numbers and coordinates, and makes sky maps
with optional annotations of objects on silicon:
```
from tessmaps import tessmaps as tm
tm.make_rect_map(sector_number, coords, names=names,
                 annotate_bools=is_transiting, title=title,
                 bkgnd_cmap='Blues', savname=savname)

```
The known exoplanets in the first science sector
[look like this](https://github.com/lgbouma/tessmaps/blob/master/results/tess_rectmap_knownplanets_sector0.png).

The clusters observed in the first science sector
[look like this](https://github.com/lgbouma/tessmaps/blob/master/results/tess_rectmap_kharchenko13_sector0.png).

You can also make all-sky maps,
[like this](https://github.com/lgbouma/tessmaps/blob/master/results/tess_pointings_radec_south_top250k.png).

Take a look at `/results/` to see similar plots for the southern hemisphere's
survey.

Finally, there's a module that converts this information to text and csv files:
```
tm.make_sector_list(sector_number, coords, names=names, savname=savname)
```
the csv files contain TIC IDs, sector properties, coordinates, and optionally
the passed names.

You can also use the shell scripts to regenerate everything. (See the docstrings).


INSTALL
----------
Currently the only install is from source. To seemlessly install the virtual
environment, use anaconda:

```
cd $SOME_DIRECTORY
git clone https://github.com/lgbouma/tessmaps
cd tessmaps
conda env create -f environment.yml -n tmaps
source activate tmaps
(tmaps) python setup.py install
```

A nasaexoplanetarchive query to crossmatch positions to planet names currently
(2018/08/06) depends on a bleeding-edge `astroquery` build. To make that:
```
cd $SOME_DIRECTORY
git clone https://github.com/astropy/astroquery
cd astroquery
source activate tmaps
(tmaps) python setup.py install
```
Then you should be ready.


IF YOU WRITE A PAPER WITH THESE TOOLS, PLEASE CITE
----------
* the TESS mission
* The TIC7.1 paper by Stassun et al for coordinates and properties of TIC objects
* Sullivan+ 2015 and Bouma+ 2017 for the focal glame geometry model
* lists appropriate to whatever you overplot:
  * if you use the knownplanets list, cite Stephen Kane
  * if you plot out the Kharchenko+ 2013 catalog clusters, cite them
* astropy
* astroquery

TODO
----------
* implement `tm.make_sector_list`
* consider making pip installable
