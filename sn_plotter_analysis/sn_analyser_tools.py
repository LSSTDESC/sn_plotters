#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 09:13:03 2024

@author: philippe.gris@clermont.in2p3.fr
"""
import glob
import pandas as pd
import h5py
from sn_tools.sn_utils import multiproc
from astropy.table import Table


def load_OS_table(dbDir, dbName, runType, season=1, fieldType='DDF'):
    """
    Function to load OS data

    Parameters
    ----------
    dbDir : str
        data location dir.
    dbName : str
        db name to process.
    runType : str
        run type (spectroz or photoz).
    season: int, optional
      season to process. The default is 1.
    fieldType : str, optional
        field type (DDF or WFD). The default is 'DDF'.

    Returns
    -------
    df : pandas df
        OS data.

    """

    fullDir = '{}/{}/{}_{}'.format(dbDir, dbName, fieldType, runType)
    search_path = '{}/SN_SN_{}_*_{}.hdf5'.format(fullDir, fieldType, season)

    fis = glob.glob(search_path)

    if len(fis) == 0:
        print('pb here files not found for path', search_path)

    assert (len(fis) > 0)

    df = pd.DataFrame()

    params = {}
    for fi in fis:
        fFile = h5py.File(fi, 'r')
        keys = list(fFile.keys())
        params['fFile'] = fFile
        dfa = multiproc(keys, params, load_os_table_multi, nproc=16)

        # idx = dfa['ebvofMW'] < 0.25
        # dfa = dfa[idx]
        df = pd.concat((df, dfa))
        # break

    return df


def load_os_table_multi(keys, params, j=0, output_q=None):
    """
    Function to load a set of SN in astropy table format

    Parameters
    ----------
    keys : list(str)
        hdf5 keys.
    params : dict
        Parameters.
    j : int, optional
        tag for multiprocessing. The default is 0.
    output_q : multiprocessing queue, optional
        Where to put the result (if not None). The default is None.

    Returns
    -------
    dict or df
        Output.

    """

    df = pd.DataFrame()

    fFile = params['fFile']
    for key in keys:
        data = Table.read(fFile, path=key)
        df = pd.concat((df, data.to_pandas()))

    if output_q is not None:
        return output_q.put({j: df})
    else:
        return df


def load_OS_df(dbDir, dbName, runType, timescale_file='year',
               timeslot=1, fieldType='DDF'):
    """
    Function to load OS data

    Parameters
    ----------
    dbDir : str
        data location dir.
    dbName : str
        db name to process.
    runType : str
       run type (spectroz or photoz).
    timescale_file : str, optional
        Time scle of the files to load. The default is 'year'.
    timeslot : list(int), optional
        Time slots to process. The default is 1.
    fieldType : str, optional
        Field type to process. The default is 'DDF'.


    Returns
    -------
    df : pandas df
        OS data.

    """

    fullDir = '{}/{}/{}_{}'.format(dbDir, dbName, fieldType, runType)
    search_path = '{}/SN_{}_*_{}_{}.hdf5'.format(
        fullDir, fieldType, timescale_file, timeslot)

    fis = glob.glob(search_path)

    if len(fis) == 0:
        print('pb: no files found for path', search_path)

    df = pd.DataFrame()

    for fi in fis:
        dfa = pd.read_hdf(fi)

        # idx = dfa['ebvofMW'] < 0.25
        # dfa = dfa[idx]

        df = pd.concat((df, dfa))
        # break
    return df


def load_DataFrame(dbDir_WFD, OS_WFD, runType='spectroz',
                   timescale_file='year', timeslots=[1], fieldType='WFD',
                   norm_factor=10):
    """
    Function to load data if pandas df

    Parameters
    ----------
    dbDir_WFD : str
        data location dir.
    OS_WFD : str
        WFD db name.
    runType : str, optional
        Run type. The default is spectroz.
    timescale_file : str, optional
       Time scale of the files to process. The default is 'year'.
    timeslots : list(int), optional
       Time slots to process. The default is [1].
    fieldType : str, optional
       Field type to process. The default is 'WFD'.

    Returns
    -------
    wfd : pandas df
        Loaded data.

    """

    wfd = pd.DataFrame()
    print('processing', OS_WFD)
    for seas in timeslots:
        wfd_seas = load_OS_df(dbDir_WFD, OS_WFD, runType=runType,
                              timescale_file=timescale_file,
                              timeslot=seas, fieldType=fieldType)
        if fieldType == 'WFD':
            wfd_seas = wfd_seas.groupby(['healpixID', timescale_file, 'field']).apply(
                lambda x: get_stat(x, norm_factor)).reset_index()
        wfd = pd.concat((wfd, wfd_seas))
        del wfd_seas

    if fieldType == 'WFD':
        print('nsn tot', wfd['nsn'].sum())
    else:
        print('nsn tot', len(wfd))

    # add a year column here
    # df_y = add_year(wfd, LSSTStart)

    return wfd


def add_year(wfd, LSSTStart):
    """
    Function to estimate the year SNe Ia have been observed

    Parameters
    ----------
    wfd : pandas df
        Data to process.
    LSSTStart : float
        LSST MJD start.

    Returns
    -------
    df_y : pandas df
        Original data + year col.

    """

    rf_phase = 35.
    df_y = pd.DataFrame()
    for y in range(1, 12):
        mjd_min = LSSTStart+(y-1)*365.
        mjd_max = LSSTStart+y*365.
        wfd['mjd_min'] = mjd_min-rf_phase*(1.+wfd['z'])
        wfd['mjd_max'] = mjd_max-rf_phase*(1.+wfd['z'])
        idx = wfd['daymax'] >= wfd['mjd_min']
        idx &= wfd['daymax'] < wfd['mjd_max']
        sel = pd.DataFrame(wfd[idx])
        sel['year'] = y
        df_y = pd.concat((df_y, sel))

    df_y = df_y.drop(columns=['mjd_min', 'mjd_max'])

    return df_y


def load_Table(dbDir_WFD, OS_WFD, runType='spectroz',
               seasons=[1], fieldType='WFD'):
    """
    Function to load data if pandas df

    Parameters
    ----------
    dbDir_WFD : str
        data location dir.
    OS_WFD : str
        WFD db name.
    runType : str, optional
        Run type. The default is spectroz.
    seasons : list(int), optional
        seasons to load. The default is [1].
    fieldType : str, optional
        Type of field to process. The default is 'WFD'.

    Returns
    -------
    wfd : pandas df
        Loaded data.

    """

    wfd = pd.DataFrame()
    print('processing', OS_WFD)
    for seas in seasons:
        wfd_seas = load_OS_table(dbDir_WFD, OS_WFD, runType=runType,
                                 season=seas, fieldType=fieldType)

        wfd = pd.concat((wfd, wfd_seas))

    # add a year column here
    # df_y = add_year(wfd, LSSTStart)

    return wfd


def get_stat(grp, norm_factor, var=['sigma_c'],
             varcut=[0.04], outvar=['nsn_cosmo']):
    """
    Function to estimate stat.

    Parameters
    ----------
    grp : pandas df
        Data to process.
    norm_factor : int
        Normalisation factor.
    var : list(str), optional
        selection vars. The default is ['sigma_c'].
    varcut : list(float), optional
        selection value. The default is [0.04].
    outvar : list(str), optional
        output variable name. The default is ['nsn_cosmo'].

    Returns
    -------
    pandas df
        Processed data.

    """

    dictout = {}

    nsn = len(grp)/norm_factor
    dictout['nsn'] = [nsn]

    idx = grp['n_epochs_bef'] >= 5
    idx &= grp['n_epochs_aft'] >= 10
    idx &= grp['n_epochs_m10_p5'] >= 5
    idx &= grp['n_epochs_phase_minus_10'] >= 2

    sel = grp[idx]

    for i, vv in enumerate(var):
        idxb = sel[vv] <= varcut[i]
        dictout[outvar[i]] = len(sel[idxb])/norm_factor

    return pd.DataFrame.from_dict(dictout)
