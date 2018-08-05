What is TESS looking at right now?
----------

This program produces a list of the top _N_ objects (stars, clusters,
galaxies) that TESS is looking at, from arbitrary input catalogs.

It tells you:
* what is being observed right now
* how many sectors it will be observed for

It also makes some sick plots to understand what's in the fields. Like this:

![sick
image](https://github.com/lgbouma/tessmaps/blob/master/results/tess_rectmap_sector0_clusters.png)

The default lists it parses are:
* Kharchenko+ 2013's cluster list.
* CTL 7.1, from filtergraph, which includes:
	* cool dwarf stars
	* hot subdwarfs
	* known planets
	* WDs

TODO:
* add Mamajek 2016's pre-Gaia association census.
* best and brightest metal-poor stars(?)
* maybe add close stars (distance less than _x_ parsecs)
* make pip installable (or otherwise save the venv)
* write install instructions
* make environment friendly

Install
----------

