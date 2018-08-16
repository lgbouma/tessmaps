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

try:
    from get_data import load_kharchenko_2013, make_prioritycut_ctl, \
            make_sublist_ctl
    from get_time_on_silicon import get_time_on_silicon
except:
    from .get_data import load_kharchenko_2013, make_prioritycut_ctl, \
            make_sublist_ctl
    from .get_time_on_silicon import get_time_on_silicon


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


def _get_TIC71_prioritycut(ctldir='/users/luke/local/TIC/CTL71/'):

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


def _get_TIC71_sublist(ctldir='/Users/luke/local/TIC/CTL71/',
                       sublist=None):
    '''
    sublist: str in ["knownplanet" , ... ]
    '''

    slpath = '../data/TIC71_{:s}cut.csv'.format(sublist)
    sltspath = '../data/TIC71_{:s}cut_tess_sectors.csv'.format(sublist)

    if not os.path.exists(slpath):
        make_sublist_ctl(datadir=ctldir, sublist=sublist)
    sldf = pd.read_csv(slpath)

    if not os.path.exists(sltspath):
        ra = np.array(sldf['RA'])
        dec = np.array(sldf['DEC'])
        coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')

        _ = get_time_on_silicon(coords)
        df = sldf.join(_, how='left')
        df.to_csv(sltspath, index=False)

    else:
        df = pd.read_csv(sltspath)

    return df


def _get_kane_knownplanets():
    '''
    sublist: str in ["knownplanet" , ... ]
    '''

    slpath = '../data/MAST_Crossmatch_CTL.csv'
    sltspath = '../data/MAST_Crossmatch_CTL_tess_sectors.csv'

    sldf = pd.read_csv(slpath)

    if not os.path.exists(sltspath):
        ra = np.array(sldf['ra'])
        dec = np.array(sldf['dec'])
        coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')

        _ = get_time_on_silicon(coords)
        df = sldf.join(_, how='left', rsuffix='kane')
        df.to_csv(sltspath, index=False)

    else:
        df = pd.read_csv(sltspath)

    return df


def get_targets(sector_number, get_kharchenko_2013=True, get_TIC71=True,
               TIC71_sublist=None, get_kane_knownplanets=None):
    '''
    For a given TESS sector, what are some good objects (stars, clusters,
    galaxies) to be aware of?

    The implemented lists to parse include:
    * Kharchenko+ 2013's cluster list.
    * CTL 7.1, from filtergraph, with sublists:
        ['bright', 'cooldwarf', 'cooldwarf;bright', 'hotsubdwarf',
         'knownplanet', 'knownplanet;bright']
    * Stephen Kane's known planet list.

    Some other ideas:
    * Mamajek 2016's pre-Gaia association census.
    * Gagne et al 2018's BANYAN association list.
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

        subcols = ['Name','d','logt','N1sr2','total_sectors_obsd']

        print('\n'+'*'*50)
        print('All Kharchenko+ 2013 clusters in sector {:d}:\n'.format(
              sector_number))

        sel = df['total_sectors_obsd']>0
        sel &= df['sector_{:d}'.format(sector_number)]>0
        print(df[subcols][sel].
              sort_values(['total_sectors_obsd','logt'],
                          ascending=[False,True]).
              to_string(index=False,col_space=12))

        csvpath = '../results/kharchenko13_sector{:d}.csv'.format(
                sector_number)
        df[subcols][sel].sort_values(
            ['total_sectors_obsd','logt'],
            ascending=[False,True]).to_csv(csvpath, index=False)
        print('\nsaved to {:s}'.format(csvpath))

        print('\n100 most-observed Kharchenko+ 2013 clusters in South:\n')

        sel = df['total_sectors_obsd']>0
        print(df[subcols][sel].
              sort_values(['total_sectors_obsd','logt'],
                          ascending=[False,True]).
              head(100).
              to_string(index=False,col_space=12))

        print('\n'+'*'*50)

    if get_TIC71:

        oksublists = ['bright', 'cooldwarf', 'cooldwarf;bright', 'hotsubdwarf',
                      'knownplanet', 'knownplanet;bright']
        if TIC71_sublist not in oksublists:
            raise AssertionError

        #typical approach: we want the thing direct from the TIC's list.
        df = _get_TIC71_sublist(sublist=TIC71_sublist)
        # np array version of e.g., ("knownplanet" in arrayname)
        sel = np.flatnonzero(np.core.defchararray.find(
                np.array(df['SPEC_LIST']).astype(str),
                TIC71_sublist)!=-1)

        #special case: use Stephen Kane's knownplanet list directly
        if get_kane_knownplanets:
            del df, sel
            df = _get_kane_knownplanets()
            sel = np.isfinite(np.array(df['ra'])) #FIXME hacky

        from plot_skymaps_of_targets import _get_knownplanet_names_transits
        names, is_transiting = _get_knownplanet_names_transits(
            df.iloc[sel], is_kane_list=get_kane_knownplanets)

        dfsel = df.iloc[sel]
        del sel
        if TIC71_sublist == 'knownplanet':
            dfsel['pl_hostname'] = names
            dfsel['is_transiting'] = is_transiting
            subcols = ['pl_hostname','is_transiting','PRIORITY',
                       'TESSMAG','RADIUS','total_sectors_obsd']
            if get_kane_knownplanets:
                subcols = ['pl_hostname','is_transiting',
                           'Tmag','total_sectors_obsd']
        else:
            dfsel = -1
            raise NotImplementedError

        print('\n'+'*'*50)
        if not get_kane_knownplanets:
            print('All {:s} from TIC7.1 in sector {:d}:\n'.
                  format(TIC71_sublist, sector_number))
        else:
            print('All {:s} from Stephen Kane\'s list in sector {:d}:\n'.
                  format(TIC71_sublist, sector_number))

        sel = dfsel['total_sectors_obsd']>0
        sel &= dfsel['sector_{:d}'.format(sector_number)]>0
        if TIC71_sublist == 'knownplanet':
            if not get_kane_knownplanets:
                print(dfsel[subcols][sel].
                      sort_values(['is_transiting','PRIORITY','TESSMAG','RADIUS'],
                                  ascending=[False,False,True,True]).
                      to_string(index=False,col_space=10))
            else:
                print(dfsel[subcols][sel].
                      sort_values(['is_transiting','Tmag'],
                                  ascending=[False,True]).
                      to_string(index=False,col_space=10))

                csvpath = '../results/kane_knownplanets_sector{:d}.csv'.format(
                        sector_number)
                dfsel[subcols][sel].sort_values(
                    ['is_transiting','Tmag'], ascending=[False,True]
                    ).to_csv(csvpath, index=False)
                print('\n saved {:s}'.format(csvpath))

        else:
            raise NotImplementedError

        if get_kane_knownplanets:
            print('\n500 most-observed {:s} from Stephen Kane\'s list in South:\n'.format(
                  TIC71_sublist))
        else:
            print('\n500 most-observed {:s} from TIC7.1 in South:\n'.format(
                  TIC71_sublist))

        sel = dfsel['total_sectors_obsd']>0
        if TIC71_sublist=='knownplanet':
            if not get_kane_knownplanets:
                print(dfsel[subcols][sel].
                      sort_values(['is_transiting','total_sectors_obsd',
                                   'PRIORITY','TESSMAG','RADIUS'],
                                  ascending=[False,False,False,True,True]).
                      head(500).
                      to_string(index=False,col_space=12))
            else:
                print(dfsel[subcols][sel].
                      sort_values(['is_transiting','total_sectors_obsd','Tmag'],
                                  ascending=[False,False,True]).
                      to_string(index=False,col_space=10))
        else:
            raise NotImplementedError

        print('\n'+'*'*50)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Generate lists of the best targets for desired catalog.')

    parser.add_argument('-sn', '--sector_number', type=int, default=None,
        help='0-12 for first year.')

    parser.add_argument('-k13', '--get_kharchenko_2013', action='store_true',
        help='Kharchenko+ 2013 cluster list.', default=False)

    parser.add_argument('-tic71', '--get_TIC71', action='store_true',
        help='Boolean to parse TIC71 special catalogs.', default=False)

    helpstr = "('bright'|'cooldwarf'|'cooldwarf;bright'|'hotsubdwarf'|"+\
              "'knownplanet'|'knownplanet;bright')"
    parser.add_argument('-ts', '--TIC71_sublist', type=str,
                        default=None, help=helpstr)

    parser.add_argument('-kane', '--get_kane_knownplanets', action='store_true',
        help='We\'re using Stephen Kane\'s knownplanet list.', default=False)

    args = parser.parse_args()

    get_targets(args.sector_number,
                get_kharchenko_2013=args.get_kharchenko_2013,
                get_TIC71=args.get_TIC71,
                TIC71_sublist=args.TIC71_sublist,
                get_kane_knownplanets=args.get_kane_knownplanets)
