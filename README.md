Overview
----------

This repo contains some tools for understanding what TESS is looking at.

One module accepts coordinates and tells you if they are observed:
```
import tessmaps as tm
df = tm.get_time_on_silicon(coords)
```
where `df` is a pandas DataFrame that tells you in which sector the observation
occurs, using the pointings given by
[tess.mit.edu/tess-pointing](https://tess.mit.edu/tess-pointing) and the focal
plane geometry used by Sullivan et al (2015) and Bouma et al (2017).

Another module accepts sector numbers and coordinates, and makes sky maps
with optional annotations of the things that are on silicon:
```
import tessmaps as tm
tm.make_rect_map(sector_number, coords, names=names,
                 annotate_bools=is_transiting, title=title,
                 bkgnd_cmap='Blues', savname=savname)

```
The known exoplanets in the first science sector look like 
[this](https://github.com/lgbouma/tessmaps/blob/master/results/tess_rectmap_knownplanets_sector0.png).
For example, the clusters observed in the first science sector look like this 
[this](https://github.com/lgbouma/tessmaps/blob/master/results/tess_rectmap_kharchenko13_sector0.png).
You can also make all-sky maps, like
[this](https://github.com/lgbouma/tessmaps/blob/master/results/tess_pointings_radec_south_top250k.png).
Plots for each sector are available [here](https://github.com/lgbouma/tessmaps/blob/master/results/).

A final module converts this information to text and csv files:
```
import tessmaps as tm
tm.make_sector_list(sector_number, coords, names=names,
                    annotate_bools=is_transiting, title=title,
                    bkgnd_cmap='Blues', savname=savname)

```
The results are output as text files
[like this](https://github.com/lgbouma/tessmaps/blob/master/results/knownplanet_sector0.txt).


Lists available by default
----------
* Kharchenko+ 2013's cluster list.
* CTL 7.1, from filtergraph, which includes:
  * known planets from Stephen Kane
	* known planets (as parsed from the SPEC\_LIST keyword. This is not perfect,
      for reasons that are being debugged.)
	* cool dwarf stars
	* hot subdwarfs
	* WDs

TODO
----------
* implement `tm.make_sector_list`
* add Gagne et al 2018's BANYAN association list.
* add Mamajek 2016's pre-Gaia association census.
* best and brightest metal-poor stars (Kevin Schlaufmann has worked on these)
* maybe add close stars (distance less than _x_ parsecs)
* make pip installable (or otherwise save the venv)

LICENSE
----------
See `LICENSE.txt`. Feel free to hack it into whatever you want. Useful pull
requests are gratefully accepted. Don't expect the exposure-time calculator to
be 100% accurate. Rely on official TESS-mission products for that.

INSTALL
----------
Current easiest install is from source, using anaconda. I'll add a `setup.py`
shortly. 

First,
```
cd $SOME_DIRECTORY
git clone https://github.com/lgbouma/tessmaps
cd tessmaps
conda env create -f environment.yml -n tmaps
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
Then you should be set.


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
