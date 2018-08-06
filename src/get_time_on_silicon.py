# -*- coding: utf-8 -*-
from __future__ import division, print_function

import numpy as np, pandas as pd

from astropy import units as u
from astropy import wcs
from astropy.coordinates import SkyCoord

def get_time_on_silicon(coords, lambda_init=315.8*u.degree, fov=24.*u.degree,
                        n_sectors=13):
    '''
    Given an array of astropy coordinates, compute
        a) which "sector numbers" the coordinates are observed in,
        b) the total time on TESS silicon each gets

    Args:
        coords: array of astropy coordinates

    Kwargs:
        lambda_init (float): the initial ecliptic longitude for TESS. Retrieved online
            2018/08/03 from tess.mit.edu/science/tess-pointing/.

        fov (float): the field of view (in degrees) for each TESS camera.

        n_segs (int): the number of "sectors" per hemisphere.

    Returns:
        a pandas DataFrame with columns:
            ra, dec, elat, elon,
            sector_1, sector_2, ..., sector_13, total_sectors_obsd

    Nomenclature:
        "Sector" means one "grouping" of 4 cameras. There are 13 sectors over
        the first TESS year. "View" means one pointing of a camera, in a
        sector. There are 26*4=104 views over the first two years.
    '''

    n_coords = len(coords)
    print('computing sector numbers for {:d} stars/objects'.format(n_coords))

    n_cameras = 4 # this will hopefully never change
    n_views = n_sectors*n_cameras

    # create dataframe that we will save boolean sector observations in.
    d = {'ra':coords.ra.value,
         'dec':coords.dec.value,
         'elon':coords.barycentrictrueecliptic.lon.value,
         'elat':coords.barycentrictrueecliptic.lat.value
        }
    df = pd.DataFrame(data=d)
    for sector in range(n_sectors):
        df['sector_{:d}'.format(sector)] = np.zeros_like(range(len(df)))

    # compute camera central coordinates for the first year, given the initial
    # longitude.
    views = []
    southern_elats = np.array([-18, -42, -66, -90])*u.degree
    for n_sector in range(n_sectors):
        this_elon = np.mod(lambda_init.to(u.deg).value +
                           n_sector*(360/n_sectors), 360)*u.degree
        for n_camera in range(n_cameras):
            this_elat = southern_elats[n_camera]

            this_coord = SkyCoord(lon=this_elon,lat=this_elat,
                                  frame='barycentrictrueecliptic')

            views.append([n_sector, n_camera, this_elon, this_elat,
                          this_coord.icrs.ra, this_coord.icrs.dec])

    views_columns = ['n_sector', 'n_camera', 'elon', 'elat', 'ra', 'dec']
    views = pd.DataFrame(views, columns=views_columns)

    # ccd info. see e.g., Huang et al, 2018. The gap size is a number inherited
    # from Josh Winn's code.
    ccd_pix = 4096*u.pix
    gap_pix = (2./0.015)*u.pix # about 133 pixels per gap
    pixel_scale = u.pixel_scale(fov/(ccd_pix + gap_pix))
    ccd_center = np.ones(2)*(ccd_pix + gap_pix) / 2
    delt = (np.ones(2)*u.pix).to(u.degree, pixel_scale)

    elon = coords.barycentrictrueecliptic.lon.value
    elat = coords.barycentrictrueecliptic.lat.value
    ra = coords.ra.value
    dec = coords.dec.value

    for ix, view in views.iterrows():

        # create new WCS object. See. e.g., Calabretta & Greisen 2002, paper
        # II, table 1. CRPIX: pixel-coordinates of reference pixel (center of
        # image). CRVAL: celestial long (RA) and lat (dec) of reference pixel.
        # CDELT: degrees per pixel at reference pixel location, i.e. the
        # coordinate scale. CTYPE: gnomic projection.
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = ccd_center.value
        w.wcs.crval = [view['elon'].value, view['elat'].value]
        ###w.wcs.crval = [view['ra'].value, view['dec'].value] # either works?
        w.wcs.cdelt = delt.value
        w.wcs.ctype = ['RA---TAN', 'DEC--TAN']

        # follow fits standards: The image pixel count starts with 1.
        x, y = w.wcs_world2pix(elon, elat, 1)
        ###x, y = w.wcs_world2pix(ra, dec, 1) # either works?

        try:
            # the extra "1" is because of 1-based image count.
            onchip = (   # lower left CCD
                         ((x > 0.0) &
                         (x < (ccd_pix.value-gap_pix.value)/2.-1.) &
                         (y > 0.0) &
                         (y < (ccd_pix.value-gap_pix.value)/2.-1.))
                     |   # lower right CCD
                         ((x > (ccd_pix.value+gap_pix.value)/2.-1.) &
                         (x < (ccd_pix.value+gap_pix.value)-1.) &
                         (y > 0.0) &
                         (y < (ccd_pix.value-gap_pix.value)/2.-1.))
                     |   # upper right CCD
                         ((x > (ccd_pix.value+gap_pix.value)/2.-1.) &
                         (x < (ccd_pix.value+gap_pix.value)-1.) &
                         (y > (ccd_pix.value+gap_pix.value)/2.-1.) &
                         (y < (ccd_pix.value+gap_pix.value)-1.))
                     |   # upper left CCD 
                         ((x > 0.0) &
                         (x < (ccd_pix.value-gap_pix.value)/2.-1.) &
                         (y > (ccd_pix.value+gap_pix.value)/2.-1.) &
                         (y < (ccd_pix.value+gap_pix.value)-1.))
                     )

        except Exception as e:

            print('failed computing onchip. error msg: {:s}'.format(e))
            return 0

        # save result to DataFrame
        df['sector_'+str(views.iloc[ix]['n_sector'])] += onchip.astype(int)

    # compute total sectors observed for each object
    sector_names = ['sector_'+str(ix) for ix in range(n_sector)]
    df['total_sectors_obsd'] = df[sector_names].sum(axis=1)

    print('computed sector numbers for {:d} stars/objects'.format(n_coords))

    return df


if __name__ == '__main__':

    np.random.seed(42)
    ras = np.random.uniform(low=0,high=359,size=420)*u.degree
    decs = np.random.uniform(low=-89,high=89,size=len(ras))*u.degree
    coords = SkyCoord(ra=ras, dec=decs)

    get_time_on_silicon(coords)
