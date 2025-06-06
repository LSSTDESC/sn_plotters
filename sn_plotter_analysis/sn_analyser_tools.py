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
import numpy as np


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
    # import time

    wfd = pd.DataFrame()
    print('processing', OS_WFD)
    for seas in timeslots:
        wfd_seas = load_OS_df(dbDir_WFD, OS_WFD, runType=runType,
                              timescale_file=timescale_file,
                              timeslot=seas, fieldType=fieldType)
        if fieldType == 'WFD':
            # time_ref = time.time()

            """
            wfd_seas = wfd_seas.groupby(['healpixID', timescale_file, 'field']).apply(
                lambda x: get_stat(x, norm_factor)).reset_index()
            """
            params = {}
            params['data'] = wfd_seas
            params['norm_factor'] = norm_factor
            params['timescale'] = timescale_file
            hpixes = wfd_seas['healpixID'].unique()
            wfd_seas = multiproc(hpixes, params, process_WFD_multi, nproc=8)

            # print('done', time.time()-time_ref)
        wfd = pd.concat((wfd, wfd_seas))
        del wfd_seas

    if fieldType == 'WFD':
        print('nsn tot', wfd['nsn'].sum())
    else:
        print('nsn tot', len(wfd))

    # add a year column here
    # df_y = add_year(wfd, LSSTStart)

    return wfd


def process_WFD_multi(hpixes, params, j, output_q=None):
    """
    multiprocessing for WFD seasons

    Parameters
    ----------
    hpixes : list(int)
        List of healpixID to process.
    params : dict
        parameter dict.
    j : int
        internal tag for multiprocessing.
    output_q : multiprocessing queue, optional
        Where to put the results. The default is None.

    Returns
    -------
    pandas df
        processed data.

    """

    data = params['data']
    norm_factor = params['norm_factor']
    timescale = params['timescale']

    idx = data['healpixID'].isin(hpixes)

    sel = data[idx]

    wfd_seas = sel.groupby(['healpixID', timescale, 'field']).apply(
        lambda x: get_stat(x, norm_factor)).reset_index()

    del sel
    del data

    if output_q is not None:
        return output_q.put({j: wfd_seas})
    else:
        return wfd_seas


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


class Estimate_NSN:
    def __init__(self, norm_factor=30,
                 rate='Hounsell', H0=70., Om=0.3,
                 minRFphaseQual=-10, maxRFphaseQual=35):
        """
        class to estimate nsn + error

        Parameters
        ----------
        norm_factor : float, optional
            Simulation normalization factor. The default is 30.
        rate : str, optional
            SN rate production. The default is 'Hounsell'.
        H0 : float, optional
            H0 parameter value. The default is 70..
        Om : float, optional
            Om parameter value. The default is 0.3.
        minRFphaseQual : float, optional
            min Rest-Frame phase quality selection. The default is -10.
        maxRFphaseQual : float, optional
            max Rest-Frame phase quality selection. The default is 35.

        Returns
        -------
        None

        """

        from sn_tools.sn_rate import SN_Rate
        self.sn_rate = SN_Rate(rate=rate,
                               H0=H0,
                               Om0=Om,
                               min_rf_phase=minRFphaseQual,
                               max_rf_phase=maxRFphaseQual)
        self.norm_factor = norm_factor

    def __call__(self, data):
        """
        Method to estimate nsn and err_nsn (using multiproc)

        Parameters
        ----------
        data : pandas df
            Data to process.

        Returns
        -------
        res : pandas df
            output data.

        """

        hpixes = data['healpixID'].unique()

        params = {}
        params['data'] = data

        from sn_tools.sn_utils import multiproc

        res = multiproc(hpixes, params, self.nsn_multiproc, 8)

        return res

    def nsn_multiproc(self, toproc, params, j=0, output_q=None):
        """
        Method to estimate nsn using multiproc

        Parameters
        ----------
        toproc : list(int)
            list of healpixIDs to process.
        params : dict
            parameters.
        j : int, optional
            Internal tag for multiprocessing. The default is 0.
        output_q : multiprocessing queue, optional
            where to put the data. The default is None.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """

        data = params['data']

        idx = data['healpixID'].isin(toproc)

        sel = data[idx]
        ccols = ['dbName', 'field', 'season', 'healpixID']
        res = sel.groupby(ccols).apply(
            lambda x: self.nsn_pixel(x), include_groups=False).reset_index()
        res['season'] = res['season'].astype(int)

        if output_q is not None:
            return output_q.put({j: res})
        else:
            return res

    def nsn_pixel(self, grp):
        """
        Method to estimate nsn per pixel/season/field/dbName

        Parameters
        ----------
        grp : pandas df
            Data to process.

        Returns
        -------
        res : pandas df
            output data.

        """

        # grab season length and survey_area
        season_length = grp['season_length'].mean()
        survey_area = grp['survey_area'].mean()

        # observed number of SN
        nsn_obs = len(grp)

        # get expected number of SN from rate
        zmin = np.min(grp['z'])
        zmax = np.max(grp['z'])
        zz, rate, err_rate, nsn, err_nsn, age_univ = self.sn_rate(
            zmin=zmin, zmax=zmax,
            duration=season_length,
            survey_area=survey_area,
            account_for_edges=True, dz=0.001)

        if len(nsn) == 0:
            res = pd.DataFrame()
        else:
            nsn_exp = int(np.cumsum(nsn)[-1]*self.norm_factor)

            if nsn_exp < 1:
                nsn_exp = 1
            # get the variance (binomial)
            p = nsn_obs/nsn_exp
            if p > 1:
                # to account for statistical fluctuations in the production
                p = 1
            var_nsn = nsn_exp*p*(1-p)

            sigma_nsn = np.sqrt(var_nsn)

            nsn_obs = nsn_obs/self.norm_factor
            err_nsn_obs = sigma_nsn/self.norm_factor

            years = grp['year'].unique()
            r = []
            for year in years:
                idx = grp['year'] == year
                sel = grp[idx]
                frac_year = len(sel)/self.norm_factor/nsn_obs
                r.append((year, nsn_obs*frac_year,
                         err_nsn_obs*frac_year))

            res = pd.DataFrame(
                r, columns=['year', 'nsn', 'err_nsn'])

        return res


def count_all(data, columns, var=['nsn'], err_var=['err_nsn']):
    """
    Function to estimate NSN and err_NSN from groupby (columns)

    Parameters
    ----------
    data : pandas df
        Data to process.
    columns : list(str)
        List of groupby columns.
    var : list(str), optional
        list of var to sum. The default is 'nsn'.
    err_var : list(str), optional
        list of var error to sum. The default is 'err_nsn'.
    Returns
    -------
    tt : pandas df
        Result.

    """

    tt = data.groupby(columns).apply(lambda x: count(
        x, var=var, err_var=err_var), include_groups=False).reset_index()

    return tt


def count(grp, var, err_var):
    """
    Function to estimate sum var, err_var

    Parameters
    ----------
    grp : pandas df
        data to process.
    var : str
        var to sum.
    err_var : str
        var error to sum.
    Returns
    -------
    res : pandas df
        Result.

    """

    dd = {}
    for vv in var:
        dd[vv] = [grp[vv].sum()]
    for err_vv in err_var:
        dd[err_vv] = [np.sqrt(grp[err_vv]**2).sum()]

    res = pd.DataFrame.from_dict(dd)

    return res


def clean_level(tt):
    """
    Function to clean the level

    Parameters
    ----------
    tt : pandas df
        Data to process.

    Returns
    -------
    tt : pandas df
        cleaned df.

    """

    tt = tt[tt.columns.drop(list(tt.filter(regex='level')))]

    return tt


def print_nsn_latex(sn_df):
    """
    Function to print latex tables of (NSN,err_NSN)

    Parameters
    ----------
    sn_df : pandas df
        Data to process.

    Returns
    -------
    None.

    """

    idx = sn_df['year'] <= 10
    sn_df = sn_df[idx]
    idxb = sn_df['year'] <= 5
    sn_df_5 = sn_df[idxb]

    res = get_nsn_pr(sn_df)
    all_years = get_nsn_pr(sn_df, cols=['dbName'])
    five_years = get_nsn_pr(sn_df_5, cols=['dbName'])

    # get the list of dbNames
    dbNames = sorted(res['dbName'].unique().tolist())
    str_db = '& '.join(dbNames)

    years = res['year'].unique()

    print('\\begin{table}[!htbp]')
    print('\\begin{center}')
    print('\caption{mycaption}\label{tab:mylabel}')
    ccols = ['|c']*len(dbNames)
    print('\\begin{tabular}{l'+''.join(ccols)+'}')
    print('\\hline')
    print('\\hline')
    linea = ' year & {} \\\\'.format(str_db)
    print(linea)
    print('\\hline')

    for year in years:
        idx = res['year'] == year
        sel = res[idx]
        mystr = '{}'.format(int(year))
        mystr += nsn_dbName(sel, dbNames)
        mystr += '\\\\'

        print(mystr)
        if year == 5:
            print('\\hline')
            mystr = '1-5'
            mystr += nsn_dbName(five_years, dbNames)
            print(mystr)
            print('\\hline')
    # 10 years
    mystr = '1-10'
    mystr += nsn_dbName(all_years, dbNames)
    print('\\hline')
    print(mystr)
    print('\\hline')
    print('\end{tabular}')
    print('\end{center}')
    print('\end{table}')


def nsn_dbName(sel, dbNames):
    """
    Function to print (NSN, err_NSN) per year

    Parameters
    ----------
    sel : pandas df
        Data to process.
    dbNames : list(str)
        list of OS to consider.

    Returns
    -------
    mystr : str
        output.

    """

    mystr = ''
    for dbName in dbNames:
        idxb = sel['dbName'] == dbName
        selb = sel[idxb]
        nsn = selb['nsn'].mean()
        err_nsn = selb['err_nsn'].mean()
        mystr += ' & {} \pm {}'.format(int(nsn), int(err_nsn))

    return mystr


def get_nsn_pr(sn_df, cols=['year', 'dbName']):
    """
    Function to grab nsn,err_nsn values

    Parameters
    ----------
    sn_df : pandas df
        Data to process.
    cols : list(str), optional
        cols to use to estimate nsn,err_nsn. The default is ['year', 'dbName'].

    Returns
    -------
    res : pandas df
        output resu.

    """

    res = count_all(sn_df, cols, var=['nsn'], err_var=['err_nsn'])

    res[['nsn', 'err_nsn']] = res[['nsn', 'err_nsn']].astype(int)

    return res
