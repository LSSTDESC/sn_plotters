#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 08:58:46 2024

@author: philippe.gris@clermont.in2p3.fr
"""
import pandas as pd
import numpy as np
from . import plt
from .sn_analyser_tools import get_stat
import os


def process_WFD_OS_nsn(conf_df, dataType, dbDir_WFD, runType,
                       timescale_file, timeslots, norm_factor,
                       outName='res_nsn_wfd.csv'):
    """
    Function to process WFD data

    Parameters
    ----------
    conf_df : pandas df
        config file.
    dataType : str
        Data type.
    dbDir_WFD : str
        Data dir.
    runType : str
        Run type.
    timescale_file : str
        Time scale (year/season)
    timeslots : list(int)
        Time slots
    norm_factor : float
        Normalization factor

    Returns
    -------
    wfd : pandas df
        Output data.

    """

    # if the results are to be processed
    if not os.path.isfile(outName):
        # get data
        nsn_wfd = get_nsn_wfd(conf_df, dataType, dbDir_WFD, runType,
                              timescale_file, timeslots, norm_factor)

        if 'level_1' in nsn_wfd.columns:
            nsn_wfd = nsn_wfd.drop(columns=['level_1'])
        nsn_wfd.to_csv(outName, index=False)

    # plot here
    plot_summary_year(outName, conf_df, timescale=timescale_file, cumul=False)
    plot_summary_year(outName, conf_df, timescale=timescale_file, cumul=True)

    plot_summary(outName, conf_df, timescale=timescale_file)
    plt.show()


def merge_ref(data, refdb='baseline_v3.4_10yrs', timescale='year'):
    """
    Function to merge data with e ref.

    Parameters
    ----------
    data : pandas df
        Data to process.
    refdb : str, optional
        Ref db Name. The default is 'baseline_v3.4_10yrs'.
    timescale : str, optional
        Timescale . The default is 'year'.

    Returns
    -------
    data : pandas df
        Data +ref.

    """

    idx = data['dbName'] == refdb
    seldb = pd.DataFrame(data[idx])

    data = data.merge(seldb, left_on=[timescale],
                      right_on=[timescale],
                      suffixes=['', '_ref'])

    return data


def plot_summary(fileName, conf_df, timescale='year', lasttime=10):
    """
    Function to plot a summary of the results (10 yrs)

    Parameters
    ----------
    fileName : str
        File with the results to process.
    conf_df : pandas df
        configuration for the plot.
    timescale : str, optional
        Timescale for the plot. The default is 'year'.
    lasttime : int, optional
        Season/year for results. The default is 10.

    Returns
    -------
    None.

    """

    data = pd.read_csv(fileName)

    refdb = 'baseline_v3.4_10yrs'

    data = merge_ref(data, refdb, timescale=timescale)

    # cumulate
    data = data.groupby(['dbName']).apply(
        lambda x: cumsumIt(x)).reset_index()

    data['delta_nsn_cosmo'] = (
        data['nsn_cosmo']-data['nsn_cosmo_ref'])/data['nsn_cosmo_ref']

    # remove refdb
    idx = data['dbName'] == refdb
    data = data[~idx]

    # choose last year
    idx = data[timescale] == lasttime

    sel = data[idx]

    fig, ax = plt.subplots(figsize=(15, 9))
    fig.subplots_adjust(bottom=0.25)
    ttit = 'ref: {} \n'.format(refdb)
    ttit += 'year {}'.format(lasttime)
    fig.suptitle(ttit)
    sel = sel.sort_values(by=['delta_nsn_cosmo'], ascending=False)
    pref = '_v3.4_10yrs'
    sel['dbName'] = sel['dbName'].str.split(pref).str[0]
    pref = 'v3.4_10yrs'
    sel['dbName'] = sel['dbName'].str.split(pref).str[0]
    ax.plot(sel['dbName'], sel['delta_nsn_cosmo'], color='k',
            linestyle='dotted', marker='o', mfc='NONE', ms=10)

    ax.grid(visible=True)
    ylab = '$\\frac{\\Delta N_{SN}}{N_{SN}}$'
    ax.set_ylabel(r'{}'.format(ylab), fontsize=30)

    plt.setp(ax.get_xticklabels(), rotation=45,
             ha="right", rotation_mode="anchor", fontsize=12)


def plot_summary_year(fileName, conf_df, timescale='year', cumul=False):
    """
    Function to plot delta nsn/nsn vs year

    Parameters
    ----------
    fileName : str
        cvs file of data to plot.
    conf_df : pandas df
        configuration for the plot.
    timescale : str, optional
        Timescale for the plot (year/season). The default is 'year'.
    cumul : bool, optional
        to plot cumul results. The default is False.

    Returns
    -------
    None.

    """

    data = pd.read_csv(fileName)

    refdb = 'baseline_v3.4_10yrs'

    data = merge_ref(data, refdb, timescale=timescale)

    if cumul:
        data = data.groupby(['dbName']).apply(
            lambda x: cumsumIt(x)).reset_index()

    data['delta_nsn_cosmo'] = (
        data['nsn_cosmo']-data['nsn_cosmo_ref'])/data['nsn_cosmo_ref']
    fig, ax = plt.subplots(figsize=(15, 9))
    fig.subplots_adjust(right=0.75)
    fig.suptitle('ref: {}'.format(refdb))
    dbNames = data['dbName'].unique().tolist()
    dbNames = sorted(dbNames)
    dbNames.remove(refdb)

    pref = 'v3.4_10yrs'

    for dbName in dbNames:
        idx = data['dbName'] == dbName
        sel = data[idx]
        idxb = conf_df['dbName_WFD'] == dbName
        selp = conf_df[idxb]
        ls = selp['ls'].values[0]
        color = selp['color'].values[0]
        marker = selp['marker'].values[0]
        dbNameb = dbName.split(pref)[0]
        if dbNameb[-1] == '_':
            dbNameb = dbNameb[:-1]
        ax.plot(sel[timescale], sel['delta_nsn_cosmo'],
                linestyle=ls, color=color, marker=marker, label=dbNameb, mfc='None')

    ax.grid(visible=True)
    ax.legend(loc='upper center',
              bbox_to_anchor=(1.20, 0.7),
              ncol=1, fontsize=12, frameon=False)
    ax.set_xlabel(r'year', fontweight='bold')
    ylab = '$\\frac{\\Delta N_{SN}}{N_{SN}}$'
    if cumul:
        ylab = '$\\frac{\\Delta \\Sigma N_{SN}}{\\Sigma N_{SN}}$'

    ax.set_ylabel(r'{}'.format(ylab), fontsize=30)


def cumsumIt(grp, timescale='year', vara='nsn_cosmo', varb='nsn_cosmo_ref'):
    """
    Function to estimate cumul results

    Parameters
    ----------
    grp : pandas df
        Data to process.
    timescale : str, optional
        Timescale (year/season). The default is 'year'.
    vara : str, optional
        var to cumul. The default is 'nsn_cosmo'.
    varb : str, optional
        var to cumul. The default is 'nsn_cosmo_ref'.

    Returns
    -------
    df : pandas df
        with cumulative results.

    """

    dictout = {}
    dictout[timescale] = grp[timescale]
    dictout[vara] = np.cumsum(grp[vara])
    dictout[varb] = np.cumsum(grp[varb])
    dictout['dbName'] = grp['dbName']
    df = pd.DataFrame.from_dict(dictout)

    return df


def get_nsn_wfd(conf_df, dataType, dbDir_WFD, runType,
                timescale_file, timeslots, norm_factor):
    """
    Function to estimate nsn from the WFD

    Parameters
    ----------
    conf_df : pandas df
        configuration parameters.
    dataType : str
        Data type to process.
    dbDir_WFD : str
        WFD dir files.
    runType : str
        run type.
    timescale_file : str
        time scale (year/season).
    timeslots : list(int)
        time slots.
    norm_factor : float
        norm factor.

    Returns
    -------
    wfd : pandas df
        Processed data.

    """

    OS_WFDs = conf_df['dbName_WFD'].unique()
    wfd = pd.DataFrame()
    from_to_load = 'from sn_plotter_analysis.sn_analyser_tools'
    mod_to_load = '{} import load_{}'.format(from_to_load, dataType)
    exec(mod_to_load)
    for OS_WFD in OS_WFDs:
        idx = conf_df['dbName_WFD'] == OS_WFD
        tt = 'load_{}(\'{}\',\'{}\',\'{}\',\'{}\',{})'.format(
            dataType, dbDir_WFD, OS_WFD, runType,
            timescale_file, timeslots)
        wfda = eval(tt)
        idx = wfda['fitstatus'] == 'fitok'
        idx &= wfda['ebvofMW'] < 0.25
        # idx &= wfda['sigma_c'] <= 0.04
        sel = wfda[idx]
        """
        tt = sel.groupby([timescale_file]).apply(
            lambda x: get_stat(x, norm_factor,
                               var=['sigma_c', 'sigma_mu'],
                               varcut=[0.04, 0.12],
                               outvar=['nsn_cosmo', 'nsn_cosmo_sigmamu'])).reset_index()
        """
        tt = sel.groupby([timescale_file]).apply(
            lambda x: get_stat(x, norm_factor,
                               var=['sigma_c'],
                               varcut=[0.04],
                               outvar=['nsn_cosmo'])).reset_index()

        tt['dbName'] = OS_WFD
        wfd = pd.concat((wfd, tt))

    return wfd
