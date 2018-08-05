# -*- coding: utf-8 -*-
from __future__ import division, print_function

import numpy as np, pandas as pd

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u

from glob import glob
import os

def make_prioritycut_ctl(datadir='/home/luke/local/TIC/CTL71/', prioritycut=0.0015):
    '''

    I downloaded the 2018/07/07 CTL direct from
        http://astro.phy.vanderbilt.edu/~oelkerrj/tic7_ctl1_20182606.tar.gz.

    It's only 2Gb, but regardless I put in on a storage drive.

    From the docs at https://filtergraph.com/tess_ctl:
        This portal was updated to reflect the CTL of TIC-7.1 on July 7, 2018.
        This Candidate Target List (CTL-7.1) is a compilation of several
        catalogs, including 2MASS, Gaia DR1, UCAC-4 & 5, Tycho-2, APASS DR9 and
        others. The CTL is the current best effort to identify stars most
        suitable for transit detection with TESS. Stars are considered for the
        CTL if they are: 1) identified as RPMJ dwarfs with greater than 2-sigma
        confidence; and 2) meet one of the following temperature/magnitude
        criteria: (TESSmag < 12 and Teff >= 5500K) or (TESSmag < 13 and Teff <
        5500K). Alternatively, a star is included in the CTL, regardless of the
        conditions above, if the star is a member of the bright star list
        (TESSmag < 6) or the specially curated cool dwarf, hot subdwarf, and
        known planet lists. Users who are interested only in the top 200K or
        400K stars may use a filter on the priority of 0.0017 and 0.0011
        respectively.  The full TIC & CTL will be available for download at
        MAST. The full machine-readable version of this CTL filtergraph portal
        is available as a comma-separated file at (above link).

    Kwargs:
        datadir, extracted should start looking like:

            luke@brik:~/local/TIC/CTL71$ tree -L 1
            .
            ├── 00-02.csv
            ├── 02-04.csv
            ├── 04-06.csv
            ├── 06-08.csv
            ├── 08-10.csv
            ├── 10-12.csv
            ├── 12-14.csv
            ├── 14-16.csv
            ├── 16-18.csv
            ├── 18-20.csv
            ├── 20-22.csv
            ├── 22-24.csv
            └── header.txt

        prioritycut: 0.0015 corresponds to top 300k or so.
    '''

    with open(datadir+'header.txt') as f:
        hdr = f.readlines()[0]
    columns = hdr.split(',')

    subcols = ['RA', 'DEC', 'TESSMAG', 'TEFF', 'PRIORITY', 'RADIUS', 'MASS',
               'CONTRATIO', 'ECLONG', 'ECLAT', 'DIST', 'TICID', 'SPEC_LIST']

    subcats = np.sort(glob(datadir+'??-??.csv'))

    for ix, subcat in enumerate(subcats):
        print(ix)
        if os.path.exists(datadir+'temp_{:d}.csv'.format(ix)):
            continue
        sc = pd.read_csv(subcat, names=columns)
        sc = sc[subcols]
        sc = sc[sc['PRIORITY']>prioritycut]

        sc.to_csv(datadir+'temp_{:d}.csv'.format(ix), index=False)

    temps = np.sort(glob(datadir+'temp_*.csv'))

    for ix, temp in enumerate(temps):
        if ix == 0:
            df = pd.read_csv(temp)
        else:
            new = pd.read_csv(temp)
            df = pd.concat([df, new])
            print('length of priorty-cut TIC list is {:d}'.format(len(df)))
        os.remove(temp)

    df.to_csv('../data/TIC71_prioritycut.csv', index=False)
    print('saved ../data/TIC71_prioritycut.csv')


def load_kharchenko_2013():
    '''
    Returns:
        close, far, df: pandas DataFrames for Kharchenko+2013 clusters, cut to
        be <500pc, <1000pc, and all the clusters.
    '''
    # Downloaded the MWSC from
    # http://cdsarc.u-strasbg.fr/viz-bin/Cat?cat=J%2FA%2BA%2F558%2FA53&target=http&
    tab = Table.read('../data/Kharchenko_2013_MWSC.vot', format='votable')

    df = tab.to_pandas()

    for colname in ['Type', 'Name', 'n_Type', 'SType']:
        df[colname] = [e.decode('utf-8') for e in list(df[colname])]

    # From Sullivan erratum:
    # For the Sun-like star, a 4 Re planet produces a transit depth of 0.13%. The
    # limiting magnitude for transits to be detectable is about I_C = 11.4 . This
    # also corresponds to K_s ~= 10.6 and a maximum distance of 290 pc, assuming no
    # extinction.
    # So, let "close" clusters be those within 500pc. "Far" clusters <1kpc.

    cinds = np.array(df['d']<500)
    close = df[cinds]
    finds = np.array(df['d']<1000)
    far = df[finds]

    N_c_r0 = int(np.sum(close['N1sr0']))
    N_c_r1 = int(np.sum(close['N1sr1']))
    N_c_r2 = int(np.sum(close['N1sr2']))
    N_f_r0 = int(np.sum(far['N1sr0']))
    N_f_r1 = int(np.sum(far['N1sr1']))
    N_f_r2 = int(np.sum(far['N1sr2']))

    type_d = {'a':'association', 'g':'globular cluster', 'm':'moving group',
              'n':'nebulosity/presence of nebulosity', 'r':'remnant cluster',
              's':'asterism', '': 'no label'}

    ntype_d = {'o':'object','c':'candidate','':'no label'}

    print('*'*50)
    print('Stats from Kharchenko+ 2013:')
    print('\nMilky Way Star Clusters (close := <500pc)'
          '\nN_clusters: {:d}'.format(len(close))+\
          '\nN_stars (in core): {:d}'.format(N_c_r0)+\
          '\nN_stars (in central part): {:d}'.format(N_c_r1)+\
          '\nN_stars (in cluster): {:d}'.format(N_c_r2))

    print('\n'+'*'*50)
    print('\nMilky Way Star Clusters (far := <1000pc)'
          '\nN_clusters: {:d}'.format(len(far))+\
          '\nN_stars (in core): {:d}'.format(N_f_r0)+\
          '\nN_stars (in central part): {:d}'.format(N_f_r1)+\
          '\nN_stars (in cluster): {:d}'.format(N_f_r2))

    print('\n'+'*'*50)
    print('\nMilky Way Star Clusters (all)'
          '\nN_clusters: {:d}'.format(len(df))+\
          '\nN_stars (in core): {:d}'.format(int(np.sum(df['N1sr0'])))+\
          '\nN_stars (in central part): {:d}'.format(int(np.sum(df['N1sr1'])))+\
          '\nN_stars (in cluster): {:d}'.format(int(np.sum(df['N1sr2']))))

    print('\n'+'*'*50)

    ####################
    # Post-processing. #
    ####################
    ra = np.array(df['RAJ2000'])
    dec = np.array(df['DEJ2000'])

    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')

    galactic_long = np.array(c.galactic.l)
    galactic_lat = np.array(c.galactic.b)
    ecliptic_long = np.array(c.barycentrictrueecliptic.lon)
    ecliptic_lat = np.array(c.barycentrictrueecliptic.lat)

    df['galactic_long'] = galactic_long
    df['galactic_lat'] = galactic_lat
    df['ecliptic_long'] = ecliptic_long
    df['ecliptic_lat'] = ecliptic_lat

    return close, far, df
