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
                        Overplot Kharchenko+ 2013 clusters on TIC7.1 rect map.
  -kp, --known_planets  Overplot Kharchenko+ 2013 clusters on TIC7.1 rect map.
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

from astropy import units as u
from astropy.coordinates import SkyCoord

##########
# helper functions
##########
def _get_knownplanet_names_transits(df):
    '''
    assign names to known planets (TIC7.1 doesn't have them). do it
    by crossmatching to NASA exoplanet archive, and taking the
    closest star as the true name. The astroquery call below, which
    is needed to know whether they transit, requires bleeding-edge
    astroquery dev branch.

    args:
        df: a dataframe of TIC71 rows and columns
    '''
    from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive
    ea = NasaExoplanetArchive.get_confirmed_planets_table(
            all_columns=True)

    ticra, ticdec = np.array(df['RA']), np.array(df['DEC'])
    neara, neadec = np.array(ea['ra']), np.array(ea['dec'])

    c_tic = SkyCoord(ra=ticra*u.deg, dec=ticdec*u.deg, frame='icrs')
    c_nea = SkyCoord(ra=neara*u.deg, dec=neadec*u.deg, frame='icrs')
    from astropy.coordinates import match_coordinates_sky
    idx_nea, d2d, _ = match_coordinates_sky(c_tic, c_nea)

    names = np.array(ea[idx_nea]['pl_hostname'])
    is_transiting = np.array(ea[idx_nea]['pl_tranflag'])
    return names, is_transiting


def _kharchenko_2013_rectmap(sector_number):

    if sector_number in [0,1,2,3,4,12]:
        overplotd = {'name':'clusters',
                     'annotateconditions':['d < 2000']}
        titlea = '\ngray: Kharchenko+13 clusters. label if d<2kpc.'
    else:
        overplotd = {'name':'clusters',
                     'annotateconditions':['d < 500', 'logt < 9']}
        titlea = '\ngray: Kharchenko+13 clusters. label if d<500pc, age<1Gyr.'

    plot_rectmaps(sector_number, seabornmap='Blues',
                  overplotd=overplotd, title_append=titlea)


def _known_planets_rectmap(sector_number):

    overplotd = {'name':'knownplanet',
                 'annotateconditions':['']}
    titlea = '\ngray: Known planets. '+\
             'Orange if observed, green if also transits.'.format(sector_number)
    plot_rectmaps(sector_number, seabornmap='Blues', overplotd=overplotd,
                  title_append=titlea)

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


def plot_yarmulkes():
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


def plot_rect(lon, lat, nsectors, size=1, title='title', savdir='../results/',
              savname=None, seabornmap='Paired', overplotd=None,
              sector_number=None, plotknownpoints=False):
    '''
    overplotd (dict): for example,
        {'name': 'clusters', 'annotateconditions': ['d < 1000', 'logt < 9']}
        The condition format must be "(key) </>/== (value)" to match the
        catalog being parsed.
    '''

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
        ksize=30
        klon = np.array([-15,15,-30,3,130])
        klat = np.array([-50, -45, -15, -80,-80])
        knsectors = np.array([5,2,3,4,12])

    # make colormap
    import seaborn as sns
    if seabornmap=='Paired':
        rgbs = sns.color_palette(seabornmap, n_colors=13, desat=0.9)
        cbarbounds = list(np.arange(0.5,14.5,1))
        bounds= list(np.arange(0.5,14.5,1))
    elif seabornmap=='Blues':
        rgbs = sns.color_palette(seabornmap, n_colors=10, desat=1.)
        rgbs = rgbs[3:]
        cbarbounds = list(np.arange(0.5,14.5,1))
        bounds= list(np.arange(0.5,4.5,1))
    else:
        raise NotImplementedError
    cmap = mpl.colors.ListedColormap(rgbs)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    if plotknownpoints:
        # sanity check: plotting the known points. are they where u expect?
        cax = ax.scatter(klon, klat, c=knsectors, s=ksize, lw=0, zorder=2,
                         cmap=cmap, norm=norm, rasterized=True,
                         transform=ccrs.PlateCarree())

    # plot the given stars that are passed (usually from TIC)
    lon[lon>180] -= 360
    cax = ax.scatter(lon, lat, c=nsectors, s=size, lw=0, zorder=2,
                     cmap=cmap, norm=norm, rasterized=True,
                     transform=ccrs.PlateCarree())
    ticlon, ticlat = lon, lat

    if type(overplotd)==dict:

        if overplotd['name']=='clusters':
            from get_targets import _get_kharchenko_2013
            df = _get_kharchenko_2013()
            lon, lat = np.array(df['ecliptic_long']), np.array(df['ecliptic_lat'])
            names = np.array(df['Name'])
        elif overplotd['name']=='knownplanet':
            from get_targets import _get_TIC71_sublist
            df = _get_TIC71_sublist(sublist=overplotd['name'])
            # np array version of e.g., ("knownplanet" in arrayname)
            sel = np.flatnonzero(np.core.defchararray.find(
                    np.array(df['SPEC_LIST']).astype(str),
                    overplotd['name'])!=-1)

            lon, lat = np.array(df['ECLONG'])[sel], np.array(df['ECLAT'])[sel]

            #FIXME
            names, is_transiting = _get_knownplanet_names_transits(df.iloc[sel])


        else:
            raise NotImplementedError

        sel = df['total_sectors_obsd']>0
        sel &= df['sector_{:d}'.format(sector_number)]>0

        # plot the positions of objects.
        lon[lon>180] -= 360
        _ = ax.scatter(lon[sel], lat[sel], c='darkorange', s=10, lw=0, zorder=4,
                       rasterized=True, transform=ccrs.PlateCarree())
        _ = ax.scatter(lon, lat, c='lightgray', s=1, lw=0, zorder=3,
                       rasterized=True, transform=ccrs.PlateCarree())

        if overplotd['name']=='knownplanet':
            sel &= is_transiting==1
            _ = ax.scatter(lon[sel], lat[sel], c='lime', s=10, lw=0,
                           zorder=5, rasterized=True,
                           transform=ccrs.PlateCarree())

        # annotate them
        subsel = sel
        conditions = overplotd['annotateconditions']
        for condition in conditions:
            if condition == '':
                continue
            cs = condition.split(' ')
            if cs[1] == '<':
                subsel &= df[cs[0]] < float(cs[2])
            elif cs[1] == '>':
                subsel &= df[cs[0]] > float(cs[2])
            elif cs[1] == '==':
                subsel &= df[cs[0]] == float(cs[2])
            else:
                raise NotImplementedError
        transform = ccrs.PlateCarree()._as_mpl_transform(ax)

        arrowprops = dict(facecolor='gray', edgecolor='gray', arrowstyle='->',
                          linewidth=0.5, connectionstyle='arc3,rad=-0.05')
        bbox = dict(facecolor='white',edgecolor='gray',
                    alpha=0.95,linewidth=0.5,pad=0.2)

        middle_elon = np.mean(ticlon[np.abs(ticlat)<45])
        diff = 45
        elon_start= np.remainder(middle_elon + diff, 360)
        elon_stop = np.remainder(middle_elon - diff, 360)
        text_elon = np.linspace(elon_start,elon_stop,len(lat[subsel]))
        text_elat = -20*np.ones_like(text_elon)
        if sector_number in [5,6,7,9,10,11]:
            text_elon = np.remainder(
                        np.linspace(elon_start,elon_stop+360,len(lat[subsel])),
                        360)

        for ix, sname, alon, alat in list(
            zip(range(len(names[subsel])), names[subsel],
                lon[subsel], lat[subsel])):

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
            fontsize='xx-small',transform=ax.transAxes,
            ha='right',va='bottom')
    ax.text(0.01,0.01,'github.com/lgbouma/tessmaps',
            fontsize='xx-small',transform=ax.transAxes,
            ha='left',va='bottom')

    ax.set_title(title, y=1.0, fontsize='x-small')

    # set circle/rect/square boundary
    import matplotlib.path as mpath
    set_circle_boundary = False
    set_rectangle_boundary = False
    if set_circle_boundary:
        phi = np.linspace(0, 2*np.pi, 100)
        center, radius = [0.5, 0.5], 0.6
        vertices = np.vstack([np.sin(phi), np.cos(phi)]).T
        circle = mpath.Path(vertices * radius + center)
        ax.set_boundary(circle, transform=ax.transAxes)
    elif set_rectangle_boundary:
        raise NotImplementedError
        # elat, elon pairs. clockwise.
        vertices = np.array([(-70,90),(-70,-150),(-5,-30),(-5,-60),(-70,90)])
        import IPython; IPython.embed()
        rect = mpath.Path(vertices)
        ax.set_boundary(rect, transform=ccrs.PlateCarree())
    else:
        pass

    fig.tight_layout()

    fig.savefig(savdir+savname.replace('pdf','png'),dpi=400,
                bbox_inches='tight')
    print('made {:s}'.format(savdir+savname.replace('pdf','png')))


def plot_rectmaps(sector_number, seabornmap='Paired', overplotd=None,
                  title_append=''):
    '''
    # sector_number: 0-12, for southern ecl. 0: July 25, 2018.
    # seabornmap: 'Paired' or 'Blues'
    # overplotd: None or Kharchenko
    overplotd (dict): for example,
        {'name': 'clusters', 'annotateconditions': ['d < 1000', 'logt < 9']}
        The condition format must be "(key) </>/== (value)" to match the
        catalog being parsed.
    '''


    df = pd.read_csv('../data/TIC71_prioritycut_tess_sectors.csv')
    ra, dec = np.array(df['RA']), np.array(df['DEC'])
    elon, elat = np.array(df['ECLONG']), np.array(df['ECLAT'])
    totsectors = np.array(df['total_sectors_obsd'])
    this_sector = np.array(df['sector_{:d}'.format(sector_number)])
    sel = this_sector > 0

    if overplotd:
        savname ='tess_rectmap_sector{:d}_{:s}.pdf'.format(
            sector_number, overplotd['name'])
    else:
        savname ='tess_rectmap_sector{:d}.pdf'.format(sector_number)

    title = 'CTL7.1 top 250k, sector {:d}.'.format(sector_number) +\
            ' tess.mit.edu/science/tess-pointing'+title_append

    plot_rect(elon[sel], elat[sel], totsectors[sel], size=0.25, title=title,
              savname=savname, seabornmap=seabornmap,
              overplotd=overplotd, sector_number=sector_number)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Make sick maps of what's visible for TESS.")

    parser.add_argument(
        '-sn', '--sector_number', type=int, default=None,
        help='0-12 for first year.')
    parser.add_argument(
        '-k13', '--kharchenko_2013', action='store_true',
        help='Overplot Kharchenko+ 2013 clusters on TIC7.1 rect map.',
        default=False)
    parser.add_argument(
        '-kp', '--known_planets', action='store_true',
        help='Overplot Kharchenko+ 2013 clusters on TIC7.1 rect map.',
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
        plot_yarmulkes()
    if args.rectmap_plain:
        assert type(args.sector_number)==int
        plot_rectmaps(args.sector_number)
    if args.kharchenko_2013:
        assert not type(args.sector_number)==int
        for sector_number in range(13):
            _kharchenko_2013_rectmap(sector_number)
    if args.known_planets:
        assert not type(args.sector_number)==int
        for sector_number in range(13):
            _known_planets_rectmap(sector_number)
