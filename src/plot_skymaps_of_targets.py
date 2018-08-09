# -*- coding: utf-8 -*-
from __future__ import division, print_function
'''
usage: plot_skymaps_of_targets.py [-h] [-sn SECTOR_NUMBER] [-k13] [-kp] [-rm]
                                  [-as] [-sc]

Make sick maps of what's visible for TESS.

optional arguments:
  -h, --help            show this help message and exit
  -sn SECTOR_NUMBER, --sector_number SECTOR_NUMBER
                        0-12 for first year.
  -k13, --kharchenko_2013
                        Overplot Kharchenko+ 2013 clusters on TESS footprint.
  -kp, --known_planets  Overplot known planets on TESS footprint.
  -rm, --rectmap_plain  Plot rect map for a given sector number, showing
                        TIC7.1 stars.
  -as, --all_sky        Plot all-sky yarmulkes.
  -sc, --sanity_check   Make a test plot to see if your coordinate system
                        works
'''
import os, argparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd, numpy as np

import tessmaps as tm

from astropy import units as u
from astropy.coordinates import SkyCoord

##########
# helper functions
##########
def _get_knownplanet_names_transits(df, is_kane_list=False):
    '''
    get names and "is_transiting" for known planets by by crossmatching
    coordinates to NASA exoplanet archive

    The astroquery call below, which is needed to know whether they transit,
    requires bleeding-edge astroquery dev branch.

    args:
        df: a dataframe with RAs and decs for which you want the corresponding
        names and "is_transiting" exoarchive column.
    '''
    from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive
    ea = NasaExoplanetArchive.get_confirmed_planets_table(
            all_columns=True)

    if not is_kane_list:
        ticra, ticdec = np.array(df['RA']), np.array(df['DEC'])
    else:
        ticra, ticdec = np.array(df['ra']), np.array(df['dec'])
    neara, neadec = np.array(ea['ra']), np.array(ea['dec'])

    c_tic = SkyCoord(ra=ticra*u.deg, dec=ticdec*u.deg, frame='icrs')
    c_nea = SkyCoord(ra=neara*u.deg, dec=neadec*u.deg, frame='icrs')
    from astropy.coordinates import match_coordinates_sky
    idx_nea, d2d, _ = match_coordinates_sky(c_tic, c_nea)

    names = np.array(ea[idx_nea]['pl_hostname'])
    is_transiting = np.array(ea[idx_nea]['pl_tranflag'])
    return names, is_transiting


def _rectmap_plain(sector_number):

    elons, elats = np.array([0,0])*u.deg, np.array([-90,-90])*u.deg
    coords = SkyCoord(lon=elons, lat=elats, frame='barycentrictrueecliptic')
    coords = coords.icrs
    savname ='tess_rectmap_sector{:d}.pdf'.format(sector_number)
    title = 'CTL7.1 top 250k, sector {:d}.'.format(sector_number) +\
            ' tess.mit.edu/science/tess-pointing'
    tm.make_rect_map(sector_number, coords, bkgnd_cmap='Paired',
                     savname=savname, title=title, ms=0)


def _kharchenko_2013_rectmap(sector_number):

    from get_targets import _get_kharchenko_2013
    df = _get_kharchenko_2013()
    ras, decs = np.array(df['ra'])*u.deg, np.array(df['dec'])*u.deg
    coords = SkyCoord(ra=ras, dec=decs, frame='icrs')
    names = np.array(df['Name'])

    if sector_number in [0,1,2,3,4,12]:
        title = 'Sector {:d}. '.format(sector_number)+\
                 'Gray: Kharchenko+13 clusters. label if d<2kpc.'
        sel = np.array(df['d']) < 2000
    else:
        title = 'Sector {:d}. '.format(sector_number)+\
                 'Gray: Kharchenko+13 clusters. label if d<500pc, age<1Gyr.'
        sel = np.array(df['d']) < 500
        sel &= np.array(df['logt']) < 9

    savname = 'tess_rectmap_kharchenko13_sector{:d}.png'.format(sector_number)

    tm.make_rect_map(sector_number, coords, names=names, annotate_bools=sel,
                     title=title, bkgnd_cmap='Blues', savname=savname)


def _known_planets_rectmap(sector_number, using_Kane_list=True):
    '''
    at the time of writing, Stephen Kane had what was probably the best list of
    known exoplanets. Cite him if you use it!
    '''

    if not using_Kane_list:
        raise NotImplementedError

    title = 'Sector {:d}. Gray: Kane\'s known planets. '.format(sector_number)+\
            'Orange: observed, green: transits.'

    savname = 'tess_rectmap_knownplanets_sector{:d}.png'.format(sector_number)

    from get_targets import _get_kane_knownplanets
    df = _get_kane_knownplanets()
    ras, decs = np.array(df['ra'])*u.deg, np.array(df['dec'])*u.deg
    coords = SkyCoord(ra=ras, dec=decs, frame='icrs')

    if using_Kane_list:
        from plot_skymaps_of_targets import _get_knownplanet_names_transits
        names, is_transiting = _get_knownplanet_names_transits(
            df, is_kane_list=True)

    tm.make_rect_map(sector_number, coords, names=names,
                     annotate_bools=is_transiting, title=title,
                     bkgnd_cmap='Blues', savname=savname)


##########
# plotting functions
##########
def plot_mwd(lon,dec,color_val,origin=0,size=3,title='Mollweide projection',
             projection='mollweide',savdir='../results/',savname='mwd_0.pdf',
             overplot_galactic_plane=True, is_tess=False, is_radec=None):

    '''
    args, kwargs:

        lon, lat are arrays of same length. they can be (RA,dec), or (ecliptic
            long, ecliptic lat). lon takes values in [0,360), lat in [-90,90],

        is_radec: mandatory. True if (RA,dec), else elong/elat.

        title is the title of the figure.

        projection is the kind of projection: 'mollweide', 'aitoff', ...

    comments: see
    http://balbuceosastropy.blogspot.com/2013/09/the-mollweide-projection.html.
    '''
    if is_radec == None:
        raise AssertionError

    # for matplotlib mollweide projection, x and y values (usually RA/dec, or
    # lat/lon) must be in -pi<x<pi, -pi/2<y<pi/2.
    # In astronomical coords, RA increases east (left on celestial charts).
    # Here, the default horizontal scale has positive to the right.

    def _shift_lon_get_x(lon, origin):
        x = np.array(np.remainder(lon+360-origin,360)) # shift lon values
        ind = x>180
        x[ind] -=360    # scale conversion to [-180, 180]
        x=-x    # reverse the scale: East to the left
        return x

    x = _shift_lon_get_x(lon, origin)

    fig = plt.figure(figsize=(6, 4.5))
    ax = fig.add_subplot(111, projection=projection, facecolor='White')

    if is_tess:
        # set up colormap
        import seaborn as sns
        rgbs = sns.color_palette('Paired', n_colors=13, desat=0.9)
        cmap = mpl.colors.ListedColormap(rgbs)
        bounds= list(np.arange(0.5,14.5,1))
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        # plot the stars
        cax = ax.scatter(np.radians(x),np.radians(dec), c=color_val, s=size,
                         lw=0, zorder=2, cmap=cmap, norm=norm, rasterized=True)

        # set up colorbar
        cbar = fig.colorbar(cax, cmap=cmap, norm=norm, boundaries=bounds,
                            fraction=0.025, pad=0.03, ticks=np.arange(13)+1,
                            orientation='vertical')
        ylabels = np.arange(1,14,1)
        cbar.ax.set_yticklabels(map(str, ylabels))
        cbar.set_label('number of pointings', rotation=270, labelpad=10)
        cbar.ax.tick_params(direction='in')

        # label each sector. this involves computing the positions first...
        views, sector_numbers, lambda_init, n_cameras = [], 13, 315.8*u.deg, 4
        southern_elats = np.array([-18, -42, -66, -90])*u.degree
        for sector_number in range(sector_numbers):
            this_elon = np.mod(lambda_init.to(u.deg).value +
                               sector_number*(360/sector_numbers), 360)*u.deg
            for n_camera in range(n_cameras):
                this_elat = southern_elats[n_camera]
                this_coord = SkyCoord(lon=this_elon,lat=this_elat,
                                      frame='barycentrictrueecliptic')
                views.append([sector_number, n_camera, this_elon, this_elat,
                              this_coord.icrs.ra, this_coord.icrs.dec])
        views_columns = ['sector_number', 'n_camera', 'elon', 'elat', 'ra', 'dec']
        views = pd.DataFrame(views, columns=views_columns)
        subsel = (views['n_camera'] == 0)
        for ix, view in views[subsel].iterrows():
            this_elon, this_elat = view['elon'], view['elat']
            sector_number = int(view['sector_number'])
            this_coord = SkyCoord(lon=this_elon,lat=this_elat,
                                  frame='barycentrictrueecliptic')
            if is_radec:
                this_x = _shift_lon_get_x(this_coord.icrs.ra.value, origin)
                this_dec = this_coord.icrs.dec.value
            else:
                this_x = _shift_lon_get_x(this_elon.value, origin)
                this_dec = this_elat.value
            ax.text(np.radians(this_x), np.radians(this_dec),
                    'S'+str(sector_number),fontsize='small', zorder=4,
                    ha='center', va='center')

    else:
        ax.scatter(np.radians(x),np.radians(dec), c=color_val, s=size,
                   zorder=2)


    if overplot_galactic_plane:

        ##########
        # make many points, and also label the galactic center. ideally you
        # will never need to follow these coordinate transformations.
        glons = np.arange(0,360)
        glats = np.zeros_like(glons)
        coords = SkyCoord(glons*u.degree, glats*u.degree, frame='galactic')
        gplane_ra, gplane_dec = coords.icrs.ra.value, coords.icrs.dec.value
        gplane_elon = coords.barycentrictrueecliptic.lon.value
        gplane_elat = coords.barycentrictrueecliptic.lat.value
        if is_radec:
            gplane_x = _shift_lon_get_x(gplane_ra, origin)
        else:
            gplane_x = _shift_lon_get_x(gplane_elon, origin)
            gplane_dec = gplane_elat
        ax.scatter(np.radians(gplane_x),np.radians(gplane_dec),
                   c='lightgray', s=2, zorder=3)
        gcenter = SkyCoord('17h45m40.04s', '-29d00m28.1s', frame='icrs')
        gcenter_ra, gcenter_dec = gcenter.icrs.ra.value, gcenter.icrs.dec.value
        gcenter_elon = gcenter.barycentrictrueecliptic.lon.value
        gcenter_elat = gcenter.barycentrictrueecliptic.lat.value
        if is_radec:
            gcenter_x = _shift_lon_get_x(np.array(gcenter_ra), origin)
        else:
            gcenter_x = _shift_lon_get_x(np.array(gcenter_elon), origin)
            gcenter_dec = gcenter_elat
        ax.scatter(np.radians(gcenter_x),np.radians(gcenter_dec),
                   c='black', s=2, zorder=4, marker='X')
        ax.text(np.radians(gcenter_x), np.radians(gcenter_dec), 'GC',
                fontsize='x-small', ha='left', va='top')
        ##########


    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder(tick_labels+360+origin,360)
    ax.set_xticklabels(tick_labels, fontsize='x-small')
    ax.set_yticklabels(np.arange(-75,75+15,15), fontsize='x-small')

    ax.set_title(title, y=1.05, fontsize='small')
    if is_radec:
        ax.set_xlabel('ra', fontsize='x-small')
        ax.set_ylabel('dec', fontsize='x-small')
    else:
        ax.set_xlabel('ecl lon', fontsize='x-small')
        ax.set_ylabel('ecl lat', fontsize='x-small')
    ax.grid(color='lightgray', linestyle='--', linewidth=0.5, zorder=-1)

    ax.text(0.99,0.01,'github.com/lgbouma/tessmaps',
            fontsize='xx-small',transform=ax.transAxes,
            ha='right',va='bottom')
    fig.tight_layout()
    #fig.savefig(savdir+savname, bbox_inches='tight')
    fig.savefig(savdir+savname.replace('pdf','png'),dpi=350,
                bbox_inches='tight')


def sanity_check():
    # sanity check: do the coords show up where they should?
    coord = np.array([(0,30), (60,-45), (240,15), (150,-75)])
    plot_mwd(coord[:,0],coord[:,1],'b', origin=0, size=4,
             title ='Test plot_mwd', savname='test_plot_mwd.pdf',
             is_radec=True)


def plot_allsky():
    '''
    make tess all-sky maps in both (ra,dec) and (elong, elat).
    verifies that the pattern look like it should.
    '''
    df = pd.read_csv('../data/TIC71_prioritycut_tess_sectors.csv')
    ra, dec = np.array(df['RA']), np.array(df['DEC'])
    nsectors = np.array(df['total_sectors_obsd'])
    sel = nsectors > 0

    plot_mwd(ra[sel], dec[sel], nsectors[sel], origin=0, size=0.1,
             title='CTL7.1 top 250k, real pointing, south, at least 1 '
             'observeration',
             savname='tess_pointings_radec_south_top250k.pdf',
             is_tess=True, is_radec=True)

    elon, elat = np.array(df['ECLONG']), np.array(df['ECLAT'])

    plot_mwd(elon[sel], elat[sel], nsectors[sel], origin=0, size=0.1,
             title='CTL7.1 top 250k, real pointing, south, at least 1 '
             'observeration',
             savname='tess_pointings_elongelat_south_top250k.pdf',
             is_tess=True, is_radec=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Make sick maps of what's visible for TESS.")

    parser.add_argument(
        '-sn', '--sector_number', type=int, default=None,
        help='0-12 for first year.')
    parser.add_argument(
        '-k13', '--kharchenko_2013', action='store_true',
        help='Overplot Kharchenko+ 2013 clusters on TESS footprint.',
        default=False)
    parser.add_argument(
        '-kp', '--known_planets', action='store_true',
        help='Overplot known planets on TESS footprint.',
        default=False)
    parser.add_argument(
        '-rm', '--rectmap_plain', action='store_true',
        help='Plot rect map for a given sector number, showing TIC7.1 stars.',
        default=False)
    parser.add_argument(
        '-as', '--all_sky', action='store_true',
        help='Plot all-sky yarmulkes.', default=False)
    parser.add_argument(
        '-sc', '--sanity_check', action='store_true',
        help='Make a test plot to see if your coordinate system works',
        default=False)

    args = parser.parse_args()

    if args.sanity_check:
        sanity_check()

    if args.all_sky:
        plot_allsky()

    if args.rectmap_plain:
        assert type(args.sector_number)==int
        _rectmap_plain(args.sector_number)

    if args.kharchenko_2013:
        assert not type(args.sector_number)==int
        for sector_number in range(13):
            _kharchenko_2013_rectmap(sector_number)

    if args.known_planets:
        assert not type(args.sector_number)==int
        for sector_number in range(13):
            _known_planets_rectmap(sector_number)
