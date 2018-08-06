What is TESS looking at right now?
----------

This program produces a list of the top _N_ objects (stars, clusters,
galaxies) that TESS is looking at, from arbitrary input catalogs.

It tells you:
* what is being observed right now
* how many sectors it will be observed for

The default lists it parses are:
* Kharchenko+ 2013's cluster list.
* CTL 7.1, from filtergraph, which includes:
	* known planets
	* cool dwarf stars
	* hot subdwarfs
	* WDs

The results can be output as text files:
[`/results/knownplanet_sector0.txt`](https://github.com/lgbouma/tessmaps/blob/master/results/knownplanet_sector0.txt)

They can also be output as plots to understand what's in the fields. For
example, the clusters observed in the first science sector look like this:

![clusters](https://github.com/lgbouma/tessmaps/blob/master/results/tess_rectmap_sector0_clusters.png)

The known exoplanets in the first science sector look like this:

![known](https://github.com/lgbouma/tessmaps/blob/master/results/tess_rectmap_sector0_knownplanet.png)

You can also make all-sky maps, like this:

![all-sky](https://github.com/lgbouma/tessmaps/blob/master/results/tess_pointings_radec_south_top250k.png)

TODO
----------
* add Gagne et al 2018's BANYAN association list.
* add Mamajek 2016's pre-Gaia association census.
* best and brightest metal-poor stars (Kevin Schlaufmann has worked on these)
* maybe add close stars (distance less than _x_ parsecs)
* make pip installable (or otherwise save the venv)
* write install instructions
* make environment friendly

LICENSE
----------
See `LICENSE.txt`. Feel free to hack it into whatever you want. Useful pull
requests are gratefully accepted. Don't expect the exposure-time calculator to
be 100% accurate. Rely on official TESS-mission products for that.
