# -*- coding: utf-8 -*-
'''
from tessmaps import make_rect_map
make_rect_map(sector_number, coords, names=names
              annotate_bools=annotate_bools)
'''
from __future__ import division, print_function

import os, argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd, numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord

from get_time_on_silicon import get_time_on_silicon

def _get_TIC_coords_count(sector_number):

    df = pd.read_csv('../data/TIC71_prioritycut_tess_sectors.csv')
    ra, dec = np.array(df['RA']), np.array(df['DEC'])
    elon, elat = np.array(df['ECLONG']), np.array(df['ECLAT'])
    totsectors = np.array(df['total_sectors_obsd'])
    this_sector = np.array(df['sector_{:d}'.format(sector_number)])
    sel = this_sector > 0

    return elon[sel], elat[sel], totsectors[sel]


def make_rect_map(sector_number, coords,
                  names=None, annotate_bools=None,
                  title=None, bkgnd_cmap='Paired',
                  savname='tess_rectmap_TEMP.png', savdir='../results/',
                  plotknownpoints=False, ms=10):
    '''
    Make a polar sky map of what TESS looks at in a given sector. The
    background is the TESS footprint for the sector; overplotted are any
    targets of interest, passed through as RAs, decs, and optionally names.

    NB. the TESS footprint, and the "on silicon" calculator differ slightly
    from Mr Tommy B's.

    args:
        sector_number: int (0-12) of TESS observing sector
        coords: np arrays of astroy SkyCoords you want to overplot against the
            TESS footprint. (decimal degrees)
    optional kwargs:
        names: names of the coordinates. especially useful if there are some
            you want annotate, using
        annotate_bools: np array with same length as elon/elat/names, with
            boolean T/F values for whether you want those elons, elats, and names
            to be annotated if they do fall on silicon. They will also get a
            highlighted color (green, instead of orange).
        plotknownpoints (bool): option if you want to see locations of
            (elon,elat) = (-15,-50),(15,-45),(-30,-15),(3,-80),(130,-80) on
            map.
        bkgnd_cmap: 'Paired' or 'Blues'.
        ms: marker size for overplotted points.
    '''

    ticelon, ticelat, tictotsectors = _get_TIC_coords_count(sector_number)

    import cartopy.crs as ccrs

    fig = plt.figure(figsize=(4,4))
    center_long = 0
    proj = ccrs.SouthPolarStereo(central_longitude=center_long,
                                 true_scale_latitude=True)
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    # the map size is set by where the scatter points are. (but try anyway).
    minlat, maxlat = -90, -10
    ax.set_extent([-180, 180, minlat, maxlat], ccrs.PlateCarree())

    if plotknownpoints:
        ksize = 30
        klon = np.array([-15,15,-30,3,130])
        klat = np.array([-50, -45, -15, -80,-80])
        knsectors = np.array([5,2,3,4,12])

    # make colormap
    import seaborn as sns
    if bkgnd_cmap=='Paired':
        rgbs = sns.color_palette(bkgnd_cmap, n_colors=13, desat=0.9)
        cbarbounds = list(np.arange(0.5,14.5,1))
        bounds= list(np.arange(0.5,14.5,1))
    elif bkgnd_cmap=='Blues':
        rgbs = sns.color_palette(bkgnd_cmap, n_colors=10, desat=1.)
        rgbs = rgbs[3:]
        cbarbounds = list(np.arange(0.5,14.5,1))
        bounds= list(np.arange(0.5,4.5,1))
    else:
        raise NotImplementedError
    cmap = mpl.colors.ListedColormap(rgbs)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    if plotknownpoints:
        cax = ax.scatter(klon, klat, c=knsectors, s=ksize, lw=0, zorder=2,
                         cmap=cmap, norm=norm, rasterized=True,
                         transform=ccrs.PlateCarree())

    # plot the TIC stars in this sector
    ticelon[ticelon>180] -= 360
    cax = ax.scatter(ticelon, ticelat, c=tictotsectors, s=0.25, lw=0, zorder=2,
                     cmap=cmap, norm=norm, rasterized=True,
                     transform=ccrs.PlateCarree())

    # get sectors where passed coords are observed
    from get_time_on_silicon import get_time_on_silicon
    df = get_time_on_silicon(coords)

    sel = df['total_sectors_obsd']>0
    sel &= df['sector_{:d}'.format(sector_number)]>0

    # plot the positions of passed coords
    elon, elat = df['elon'], df['elat']
    elon[elon>180] -= 360
    _ = ax.scatter(elon[sel], elat[sel], c='darkorange', s=ms, lw=0, zorder=4,
                   rasterized=True, transform=ccrs.PlateCarree())
    _ = ax.scatter(elon, elat, c='lightgray', s=1, lw=0, zorder=3,
                   rasterized=True, transform=ccrs.PlateCarree())

    if type(annotate_bools)==np.ndarray:

        # first, highlight annotated coords
        sel &= annotate_bools
        _ = ax.scatter(elon[sel], elat[sel], c='lime', s=ms, lw=0,
                       zorder=5, rasterized=True,
                       transform=ccrs.PlateCarree())

        subsel = sel
        # calculate the positions of annotation labels
        middle_elon = np.mean(ticelon[np.abs(ticelat)<45])
        diff = 45
        elon_start= np.remainder(middle_elon + diff, 360)
        elon_stop = np.remainder(middle_elon - diff, 360)
        text_elon = np.linspace(elon_start,elon_stop,len(elat[subsel]))
        text_elat = -20*np.ones_like(text_elon)
        if sector_number in [4,5,6,7,9,10,11,12]:
            text_elon = np.remainder(
                        np.linspace(elon_start,elon_stop+360,len(elat[subsel])),
                        360)

        transform = ccrs.PlateCarree()._as_mpl_transform(ax)
        arrowprops = dict(facecolor='gray', edgecolor='gray', arrowstyle='->',
                          linewidth=0.5, connectionstyle='arc3,rad=-0.05')
        bbox = dict(facecolor='white',edgecolor='gray',
                    alpha=0.95,linewidth=0.5,pad=0.2)

        for ix, sname, alon, alat in list(
            zip(range(len(names[subsel])), names[subsel],
                elon[subsel], elat[subsel])):

            ax.annotate(sname, xy=(alon,alat), xycoords=transform,
                        xytext=(text_elon[ix], text_elat[ix]),
                        textcoords=transform, ha='center', va='top',
                        arrowprops=arrowprops, bbox=bbox, fontsize='xx-small',
                        zorder=4)
    else:
        pass

    # set up colorbar
    cbar = fig.colorbar(cax, cmap=cmap, norm=norm, boundaries=cbarbounds,
                        fraction=0.04, pad=0.03, ticks=np.arange(13)+1,
                        orientation='vertical')
    ylabels = np.arange(1,14,1)
    cbar.ax.set_yticklabels(map(str, ylabels))
    cbar.set_label('number of pointings', rotation=270, labelpad=10)
    cbar.ax.tick_params(direction='in')

    # make grid lines and label the spiral
    ax.gridlines(linewidth=0.5, linestyle='--', color='lightgray', zorder=-1)
    lonlabels = np.arange(-120,180+60,60)
    latlabels = -10*np.ones_like(lonlabels)
    ix = 0
    for lon, lat in list(zip(lonlabels, latlabels)):
        if lon >= 0:
            ax.text(lon, lat-ix*10, str('({:d},{:d})'.format(lon,lat-ix*10)),
                    fontsize='xx-small', transform=ccrs.PlateCarree(),
                    ha='center', va='center')
        else:
            ax.text(lon, lat-ix*10,
                    str('({:d},{:d})'.format(lon+360,lat-ix*10)),
                    fontsize='xx-small', transform=ccrs.PlateCarree(),
                    ha='center', va='center')
        ix += 1
    ax.text(0.99,0.01,'ecliptic coords',
            fontsize='xx-small',transform=ax.transAxes, ha='right',va='bottom')
    ax.text(0.01,0.01,'github.com/lgbouma/tessmaps',
            fontsize='xx-small',transform=ax.transAxes, ha='left',va='bottom')

    ax.set_title(title, y=1.0, fontsize='xx-small')

    fig.tight_layout()

    fig.savefig(savdir+savname.replace('pdf','png'),dpi=400,
                bbox_inches='tight')
    print('made {:s}'.format(savdir+savname.replace('pdf','png')))
