# -*- coding: utf-8 -*-
from __future__ import division, print_function

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd, numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord

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
        views, n_sectors, lambda_init, n_cameras = [], 13, 315.8*u.deg, 4
        southern_elats = np.array([-18, -42, -66, -90])*u.degree
        for n_sector in range(n_sectors):
            this_elon = np.mod(lambda_init.to(u.deg).value +
                               n_sector*(360/n_sectors), 360)*u.deg
            for n_camera in range(n_cameras):
                this_elat = southern_elats[n_camera]
                this_coord = SkyCoord(lon=this_elon,lat=this_elat,
                                      frame='barycentrictrueecliptic')
                views.append([n_sector, n_camera, this_elon, this_elat,
                              this_coord.icrs.ra, this_coord.icrs.dec])
        views_columns = ['n_sector', 'n_camera', 'elon', 'elat', 'ra', 'dec']
        views = pd.DataFrame(views, columns=views_columns)
        subsel = (views['n_camera'] == 0)
        for ix, view in views[subsel].iterrows():
            this_elon, this_elat = view['elon'], view['elat']
            n_sector = int(view['n_sector'])
            this_coord = SkyCoord(lon=this_elon,lat=this_elat,
                                  frame='barycentrictrueecliptic')
            if is_radec:
                this_x = _shift_lon_get_x(this_coord.icrs.ra.value, origin)
                this_dec = this_coord.icrs.dec.value
            else:
                this_x = _shift_lon_get_x(this_elon.value, origin)
                this_dec = this_elat.value
            ax.text(np.radians(this_x), np.radians(this_dec),
                    'S'+str(n_sector),fontsize='small', zorder=4,
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
              savname='temp.pdf', seabornmap='Paired',
              overplottype=None, sector_number=None):

    import cartopy.crs as ccrs

    fig = plt.figure(figsize=(4,4))
    center_long = 0
    proj = ccrs.SouthPolarStereo(central_longitude=center_long,
                                 true_scale_latitude=True)
    ax = fig.add_subplot(1, 1, 1, projection=proj)
    minlat, maxlat = -90, -10
    # minlat actually does nothing. the map size is set by where the scatter
    # points are.
    ax.set_extent([-180, 180, minlat, maxlat], ccrs.PlateCarree())

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

    ## sanity check: plot the known ones
    #cax = ax.scatter(klon, klat, c=knsectors, s=ksize, lw=0, zorder=2,
    #                 cmap=cmap, norm=norm, rasterized=True,
    #                 transform=ccrs.PlateCarree())

    # plot the stars
    lon[lon>180] -= 360
    cax = ax.scatter(lon, lat, c=nsectors, s=size, lw=0, zorder=2,
                     cmap=cmap, norm=norm, rasterized=True,
                     transform=ccrs.PlateCarree())

    if overplottype=='clusters':
        from get_targets import _get_kharchenko_2013
        df = _get_kharchenko_2013()

        subcols = ['Name','d','logt','N1sr2','total_sectors_obsd']
        lon, lat = np.array(df['ecliptic_long']), np.array(df['ecliptic_lat'])
        names = np.array(df['Name'])

        sel = df['total_sectors_obsd']>0
        sel &= df['sector_{:d}'.format(sector_number)]>0

        print('\nin sector {:d} observing:\n'.format(sector_number))
        print(df[subcols][sel].
              sort_values(['total_sectors_obsd','logt'],ascending=[False,True]).
              to_string(index=False,col_space=12))
        print('\n')

        # plot the points
        lon[lon>180] -= 360
        _ = ax.scatter(lon[sel], lat[sel], c='darkorange', s=10, lw=0, zorder=4,
                       rasterized=True, transform=ccrs.PlateCarree())
        _ = ax.scatter(lon, lat, c='lightgray', s=1, lw=0, zorder=3,
                       rasterized=True, transform=ccrs.PlateCarree())

        # annotate them
        subsel = sel
        subsel &= df['d'] < 2000     # within 2kpc
        transform = ccrs.PlateCarree()._as_mpl_transform(ax)

        arrowprops = dict(facecolor='gray', edgecolor='gray', arrowstyle='->',
                          linewidth=0.5, connectionstyle='arc3,rad=-0.05')
        bbox = dict(facecolor='white',edgecolor='gray',
                    alpha=0.95,linewidth=0.5,pad=0.2)

        text_elon = np.linspace(0,270,len(lat[subsel]))
        text_elat = -20*np.ones_like(text_elon)

        for ix, sname, alon, alat in list(
            zip(range(len(names[subsel])), names[subsel],
                lon[subsel], lat[subsel])):

            ax.annotate(sname, xy=(alon,alat), xycoords=transform,
                        xytext=(text_elon[ix], text_elat[ix]),
                        textcoords=transform, ha='center', va='top',
                        arrowprops=arrowprops, bbox=bbox, fontsize='xx-small',
                        zorder=4)

        #ax.annotate('test', xy=(0,-30), xycoords=transform, xytext=(0,-10),
        #            textcoords=transform, ha='center', va='top',
        #            arrowprops=arrowprops, bbox=bbox, fontsize='xx-small')

    elif overplottype:
        raise NotImplementedError

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


def plot_rectmaps(n_sector, seabornmap='Paired', overplottype=None,
                  sector_number=None, title_append=''):
    # n_sector: e.g., 0, if it's starting July 25, 2018.
    # seabornmap: 'Paired' or 'Blues'
    # overplottype: None or Kharchenko
    # sector_number: 0-12, for southern ecl

    df = pd.read_csv('../data/TIC71_prioritycut_tess_sectors.csv')
    ra, dec = np.array(df['RA']), np.array(df['DEC'])
    elon, elat = np.array(df['ECLONG']), np.array(df['ECLAT'])
    totsectors = np.array(df['total_sectors_obsd'])
    this_sector = np.array(df['sector_{:d}'.format(n_sector)])
    sel = this_sector > 0

    if overplottype:
        savname ='tess_rectmap_sector{:d}_{:s}.pdf'.format(n_sector, overplottype)
    else:
        savname ='tess_rectmap_sector{:d}.pdf'.format(n_sector)

    title = 'CTL7.1 top 250k, sector {:d}.'.format(n_sector) +\
            ' tess.mit.edu/science/tess-pointing'+title_append

    plot_rect(elon[sel], elat[sel], totsectors[sel], size=0.25, title=title,
              savname=savname, seabornmap=seabornmap,
              overplottype=overplottype, sector_number=sector_number)


if __name__ == '__main__':

    n_sector=0

    #sanity_check()
    plot_yarmulkes()

    plot_rectmaps(n_sector)

    plot_rectmaps(
        n_sector, seabornmap='Blues', overplottype='clusters', sector_number=0,
        title_append='\ngray: Kharchenko+13 clusters. labelled if d<2kpc.')
