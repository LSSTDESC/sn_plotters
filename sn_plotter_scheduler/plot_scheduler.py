#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 09:33:06 2025

@author: philippe.gris@clermont.in2p3.fr
"""
from . import plt
import pandas as pd
import numpy as np
from datetime import datetime
from astropy.time import Time
from sn_tools.sn_obs import season

plt.rcParams["axes.labelsize"] = "large"
plt.rcParams["axes.linewidth"] = 2.0
plt.rcParams["xtick.major.size"] = 8
plt.rcParams["ytick.major.size"] = 8
plt.rcParams["ytick.minor.size"] = 5
plt.rcParams["xtick.labelsize"] = "large"
plt.rcParams["ytick.labelsize"] = "large"

plt.rcParams["figure.figsize"] = (12, 8)
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.titleweight'] = 'bold'
# plt.rcParams['axes.facecolor'] = 'blue'
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'
# the line width around the marker symbol
plt.rcParams['lines.markeredgewidth'] = 0.3
plt.rcParams['lines.markersize'] = 5  # markersize, in points
plt.rcParams['grid.alpha'] = 0.75  # transparency, between 0.0 and 1.0
plt.rcParams['grid.linestyle'] = '-'  # simple line
plt.rcParams['grid.linewidth'] = 0.4  # in points
plt.rcParams['font.size'] = 13


def add_info(df):
    """
    Function to add infos about the df

    Parameters
    ----------
    df : pandas df
        Data to process.

    Returns
    -------
    df : pandas df
        Output data.

    """

    # add the date
    dfb = pd.DataFrame({'year': df['year'].to_list(),
                        'month': df['month'].to_list(),
                        'day': df['day'].to_list()})

    df['date'] = pd.to_datetime(dfb).to_list()

    # estimate Nvisits corresponding to obs_time
    df['nvisits'] = (df['obs_duration [h]']*3600/30).astype(int)

    # add seasons
    obs = df.to_records(index=False)

    targets = np.unique(obs['target'])

    res = None
    for tt in targets:
        idx = obs['target'] == tt
        sel = obs[idx]
        sel_seas = season(sel, season_gap=20, mjdCol='mjd')
        if res is None:
            res = sel_seas
        else:
            res = np.concatenate((res, sel_seas))

    df = pd.DataFrame.from_records(res)
    df = df.rename(columns={'period': 'season'})

    mjd_min = df['mjd'].min()
    df['night'] = (df['mjd']-mjd_min+1).astype(int)

    return df


def plot(df, fields=['COSMOS']):
    """
    Function to plot observing time vs time

    Parameters
    ----------
    df : pandas df
        Data to process.
    fields: list(str)
        List of fields to display.

    Returns
    -------
    None.

    """

    mjd_min = df['mjd'].min()
    ymax = df['nvisits'].max()+10

    for tt in fields:
        fig, ax = plt.subplots(figsize=(14, 8))
        fig.suptitle(tt)
        idx = df['target'] == tt
        sel = pd.DataFrame(df[idx])
        """
        ax.plot(sel['date'], sel['obs_duration [h]'],
                marker='.', linestyle='None')
        """
        ax.plot(sel['date'], sel['nvisits'],
                marker='.', linestyle='None', color='b')

        seasons = sel['season'].unique()

        for seas in seasons:
            ii = sel['season'] == seas
            selb = sel[ii]
            ax.plot(selb['date'], selb['nvisits'],
                    marker='.', linestyle='None')

    for nn in range(11):
        xv = mjd_min+(nn*365)
        ti = Time('{}'.format(xv), format='mjd')
        print('bbb', ti.ymdhms)
        ddate = datetime(ti.ymdhms[0], ti.ymdhms[1], ti.ymdhms[2])
        ax.plot([ddate]*2, [0, ymax], color='r', lw=2)

    ymin = 0.

    ax.set_ylim([ymin, ymax])
    ax.grid(visible=True)
    ax.set_xlabel(r'Time [year]')
    ax.set_ylabel(r'$N_{visits}$')


def ana_season(sel_df):
    """
    Function to analyze a season

    Parameters
    ----------
    sel_df : pandas df
        Data to process.

    Returns
    -------
    pandas df
        output data.

    """
    tti = np.arange(0., 4.5, 0.02)

    var = 'obs_duration [h]'
    r = []
    for i in range(len(tti)-1):
        io = sel_df[var] >= tti[i]
        io &= sel_df[var] < tti[i+1]
        ssel = sel_df[io]
        sl = ssel['mjd'].max()-ssel['mjd'].min()
        r.append((np.mean([tti[i], tti[i+1]]), sl))

    res = np.rec.fromrecords(r, names=['observable_time [h]', 'season_length'])

    return pd.DataFrame.from_records(res)


def smooth_It(vals, xvar, yvar, kk=3):
    """
    Function to smooth a set of data

    Parameters
    ----------
    vals : array
        Data to process.
    xvar : str
        x-axis variable.
    yvar : str
        y-axis variable.
    kk : int, optional
        smoothing parameter. The default is 3.

    Returns
    -------
    xnew : array
        xnew.
    spl_smooth : array
        y smoothed values.

    """

    from scipy.interpolate import make_interp_spline, UnivariateSpline
    xmin, xmax = np.min(vals[xvar]), np.max(vals[xvar])
    xnew = np.linspace(xmin, xmax, 100)
    spl = make_interp_spline(
        vals[xvar], vals[yvar], k=kk)  # type: BSpline
    spl = UnivariateSpline(vals[xvar], vals[yvar], k=kk)
    spl.set_smoothing_factor(0.5)
    spl_smooth = spl(xnew)
    return xnew, spl_smooth


def fa(x):
    return x*3600./30.


def fb(x):

    return x*3600./30.


def plot_season_length(res_season):
    """
    Function to plot season length vs obs time/Nvisits

    Parameters
    ----------
    res_season : pandas df
        Data to plot.

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots()
    # axb = ax.twiny()
    dfb = pd.DataFrame()
    tta = ['COSMOS', 'XMM_LSS', 'ELAISS1', 'ECDFS', 'EDFS_a']
    lcol = ['g', 'm', 'r', 'orange', 'b']
    lst = ['solid', 'dotted', 'dashdot', 'dashed', 'dashdot']
    labs = ['COSMOS', 'XMM-LSS', 'ELAIS', 'CDFS', 'EDFS$_{a,b}$']

    dcol = dict(zip(tta, lcol))
    dlst = dict(zip(tta, lst))
    dlabs = dict(zip(tta, labs))

    for tt in tta:
        idx = res_season['target'] == tt
        idx &= res_season['season_length'] > 50
        sel_seas = res_season[idx]
        # sel_seas['nvisits'] = (sel_seas['observable_time [h]']*3600/30).astype(int)
        dfb = pd.concat((dfb, sel_seas))

        ax.plot(sel_seas['observable_time [h]'],
                sel_seas['season_length'], marker='None',
                label=dlabs[tt], color=dcol[tt], linestyle=dlst[tt], lw=2)
        # axb.plot(sel_seas['nvisits'],
        #         sel_seas['season_length'])
        # x, y = smooth_It(sel_seas, 'observable_time [h]', 'season_length')
        # ax.plot(x, y, marker='.', label=tt)

    ax.grid(visible=True)

    ymin = dfb['season_length'].min()-5.
    ymax = dfb['season_length'].max()+5.
    ax.plot([0, 4], [180.]*2, color='k', lw=2)
    ax.set_xlim([0, 4])
    ax.set_ylim([ymin, ymax])
    ax.legend()
    ax.set_xlabel(r'Observable time [h]')
    ax.set_ylabel(r'Max season length [days]')
    ax.text(0.5, 185., 'season length = 6 months')
    secax = ax.secondary_xaxis('top', functions=(fa, fb))
    secax.set_xlabel('N$_{visits}$')
