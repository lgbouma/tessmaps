# -*- coding: utf-8 -*-
'''
usage: get_targets.py [-h] [-sn SECTOR_NUMBER] [-k13] [-tic71]

Generate lists of the best targets for desired catalog.

optional arguments:
  -h, --help            show this help message and exit
  -sn SECTOR_NUMBER, --sector_number SECTOR_NUMBER
                        0-12 for first year.
  -k13, --get_kharchenko_2013
                        Kharchenko+ 2013 cluster list.
  -tic71, --get_TIC71   TIC71 special catalog lists.
'''
from __future__ import division, print_function

import os, argparse
import numpy as np, pandas as pd

from astropy import units as u
from astropy.coordinates import SkyCoord

from get_data import load_kharchenko_2013, make_prioritycut_ctl
from get_time_on_silicon import get_time_on_silicon

def _get_kharchenko_2013():

    k2013path = '../data/Kharchenko_2013_MWSC_tess_sectors.csv'

    if not os.path.exists(k2013path):
        _, _, k_all = load_kharchenko_2013()

        ra = np.array(k_all['RAJ2000'])
        dec = np.array(k_all['DEJ2000'])
        coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')

        _ = get_time_on_silicon(coords)

        k_all = k_all.reset_index()
        df = k_all.join(_, how='left')

        df.to_csv(k2013path, index=False)

    else:
        df = pd.read_csv(k2013path)

    return df


def _get_TIC71(ctldir='/home/luke/local/TIC/CTL71/'):

    pcpath = '../data/TIC71_prioritycut.csv'
    pctspath = '../data/TIC71_prioritycut_tess_sectors.csv'

    if not os.path.exists(pcpath):
        make_prioritycut_ctl(datadir=ctldir)
    pc_ctl = pd.read_csv(pcpath)

    if not os.path.exists(pctspath):

        ra = np.array(pc_ctl['RA'])
        dec = np.array(pc_ctl['DEC'])
        coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')

        _ = get_time_on_silicon(coords)

        df = pc_ctl.join(_, how='left')

        df.to_csv(pctspath, index=False)

    else:
        df = pd.read_csv(pctspath)

    return df


def get_targets(sector_number, get_kharchenko_2013=True, get_TIC71=True):
    '''
    For a given TESS sector, what are some good objects (stars, clusters,
    galaxies) to be aware of?

    The default lists it parses are:
    * Kharchenko+ 2013's cluster list.
    * CTL 7.1, from filtergraph, which includes:
        cool dwarf stars. bright stars. hot subdwarfs. known planets.
        "knownplanet;bright" (might want to read the docs).

    Some other ideas:
    * Mamajek 2016's pre-Gaia association census.
    * WDs
    * best and brightest metal-poor stars(?)
    * close stars (distance less than _x_ parsecs)

    ----------
    Args:
        sector_number (int): indexed starting with ZERO because we are not
        anarchists.

    '''

    assert sector_number >= 0
    if sector_number > 12:
        raise NotImplementedError

    if get_kharchenko_2013:
        df = _get_kharchenko_2013()
        print('\n'+'*'*50)

        print('All Kharchenko+ 2013 clusters in sector {:d}:\n'.format(
              sector_number))

        subcols = ['Name','d','logt','N1sr2','total_sectors_obsd']

        sel = df['total_sectors_obsd']>0
        sel &= df['sector_{:d}'.format(sector_number)]>0

        print(df[subcols][sel].
              sort_values(['total_sectors_obsd','logt'],ascending=[False,True]).
              to_string(index=False,col_space=12))

        print('\n100 most-observed Kharchenko+ 2013 clusters in South:\n')

        sel = df['total_sectors_obsd']>0
        print(df[subcols][sel].
              sort_values(['total_sectors_obsd','logt'],ascending=[False,True]).
              head(100).
              to_string(index=False,col_space=12))

        print('\n'+'*'*50)

    if get_TIC71:
        df = _get_TIC71()

        speclists = np.unique(df[~df['SPEC_LIST'].isnull()]['SPEC_LIST'])

        for speclist in speclists:
            print('\n'+'*'*50)
            #FIXME here. 
            print


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Generate lists of the best targets for desired catalog.')

    parser.add_argument('-sn', '--sector_number', type=int, default=None,
        help='0-12 for first year.')

    parser.add_argument('-k13', '--get_kharchenko_2013', action='store_true',
        help='Kharchenko+ 2013 cluster list.', default=False)

    parser.add_argument('-tic71', '--get_TIC71', action='store_true',
        help='TIC71 special catalog lists.', default=False)

    args = parser.parse_args()

    if args.get_TIC71:
        raise NotImplementedError

    get_targets(args.sector_number,
                get_kharchenko_2013=args.get_kharchenko_2013,
                get_TIC71=args.get_TIC71)
