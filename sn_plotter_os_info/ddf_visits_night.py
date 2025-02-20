#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 10:01:53 2025

@author: philippe.gris@clermont.in2p3.fr
"""
import numpy as np
from . import plt
import pandas as pd
from sn_analysis.sn_calc_plot import bin_it


def analyze_simu_exp(data):
    """
    Function to analyze sim+expected data

    Parameters
    ----------
    data : pandas df
        Data to process.

    Returns
    -------
    None.

    """

    # add new columns as diff

    for vv in ['nvisits', 'u', 'g', 'r', 'i', 'z', 'y']:
        data['diff_{}'.format(vv)] = data['{}_exp'.format(vv)] - data[vv]
        data['ratio_{}'.format(vv)] = data['{}_exp'.format(vv)] / data[vv]

    print(data['ratio_u'])
    data = data.replace([np.inf, -np.inf], 0)
    print(data['ratio_u'])

    # plot stat results
    # ana_plot_stat(data)

    # analysis of cases when thee number of visits exceeds expectation
    plot_stat_visits_vs_exp(data)

    # plot obs time vs night
    # plot_obs_time_night(data)


def plot_obs_time_night(data):
    """
    Function to plot obs time and nddf per night of observation

    Parameters
    ----------
    data : pandas df
        Data to process.

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots(nrows=2, figsize=(12, 8))
    fig.subplots_adjust(hspace=0)
    print(data.columns)
    fig.suptitle(data['dbName'].unique()[0])
    obs_time = data.groupby(['night']).apply(
        lambda x: pd.DataFrame(
            {'obs_time [h]': [x['nvisits'].sum()*30./3600.],
             'nddf': [len(x['target_name'].unique())]
             })).reset_index()

    ax[0].plot(obs_time['night'], obs_time['nddf'], 'k.')

    ax[0].set_xticklabels([])
    ax[0].set_ylabel(r'N$_{DDF}$')

    ax[1].plot(obs_time['night'], obs_time['obs_time [h]'], 'k.')
    print(obs_time)

    ax[1].set_xlabel(r'night')
    ax[1].set_ylabel(r'DDF obs. time [h]')

    night_max = obs_time['night'].max()
    tdays = np.arange(1, night_max, 365.)

    colors = ['blue', 'orange', 'green', 'crimson', 'lightsalmon',
              'grey', 'yellow', 'magenta', 'purple', 'cyan']
    print(len(tdays))

    for i in range(2):
        ax[i].grid(visible=True)
        ax[i].set_xlim([0, night_max])
        ax[i].set_ylim([0, None])
        # ax[i].axvspan(1, 365, facecolor='red', alpha=0.25)
        for j in range(len(tdays)):
            ax[i].axvspan(tdays[j], tdays[j]+365.,
                          facecolor=colors[j], alpha=0.25)

    plt.show()


def plot_stat_visits_vs_exp(data):

    idx = data['diff_nvisits'] < 0.
    idx &= data['target_name'] == 'DD:COSMOS'
    sel = data[idx]

    print(sel)

    dd = {}
    for b in 'ugrizy':
        rb = bin_it(sel, 'ratio_{}'.format(b), bins=np.arange(
            0, 1.1, 0.1), norm_factor=1, outvar='frac')
        rb['frac'] /= rb['frac'].sum()
        dd[b] = rb

    fig, ax = plt.subplots(figsize=(12, 8))
    fig.suptitle(sel['dbName'].unique()[0])

    for key, vals in dd.items():
        ax.plot(vals['ratio_{}'.format(key)], 100.*vals['frac'])

    ax.grid(visible=True)
    ax.set_ylabel(r'Fraction of night [%]')
    ax.set_xlabel(r'$\frac{N_{visits}^{exp}}{N_{visits}^{obs}}$')
    ax.set_xlim([0, None])
    ax.set_ylim([0, None])
    plt.show()

    ax.hist(sel['ratio_nvisits'], histtype='step', bins=20)

    for b in 'ugrizy':
        fig, ax = plt.subplots()
        fig.suptitle('{}-band'.format(b))
        vvar = 'ratio_{}'.format(b)
        ax.hist(sel[vvar], histtype='step', bins=20)

    plt.show()


def ana_plot_stat(data):
    """
    Function to make some stat and to make plots

    Parameters
    ----------
    data : pandas df
        Data to process.

    Returns
    -------
    None.

    """

    res_stat = data.groupby(['target_name', 'season', 'DD_type']).apply(
        lambda x: stat_simu_exp(x)).reset_index()

    print(res_stat)

    plot_stat(res_stat)


def stat_simu_exp(grp):
    """
    Function to get some stat on simu+expected os data

    Parameters
    ----------
    grp : pandas df
        Data to process.

    Returns
    -------
    rr : pandas df
        output data.

    """

    dict_frac = {}
    nnights = len(grp)
    for vv in ['nvisits', 'u', 'g', 'r', 'i', 'z', 'y']:
        myvar = 'diff_{}'.format(vv)
        idxa = grp[myvar] >= 1
        idxb = grp[myvar] <= -1
        idxc = np.abs(grp[myvar]) < 0.5
        dict_frac['{}_missing'.format(vv)] = [len(grp[idxa])/nnights]
        dict_frac['{}_excess'.format(vv)] = [len(grp[idxb])/nnights]
        dict_frac['{}_perfect'.format(vv)] = [len(grp[idxc])/nnights]

    rr = pd.DataFrame.from_dict(dict_frac)

    return rr


def plot_stat(data, prefix='DD:'):
    """
    Function to plot stat results

    Parameters
    ----------
    data : pandas df
        Data to process.
    prefix : str, optional
        prefix for DD fieldnames. The default is 'DD:'.

    Returns
    -------
    None.

    """

    fields = data['target_name'].unique()
    categ = ['perfect', 'missing', 'excess']
    value = ['=', '<', '>']
    fields = ['COSMOS', 'XMM_LSS', 'ELAISS1', 'ECDFS', 'EDFS_a', 'EDFS_b']
    fields = list(map(lambda x: prefix + x, fields))
    colors = ['r', 'b', 'k', 'orange', 'm', 'g']
    marks = ['o', 's', '*', '^', 'v', '>']
    dict_col = dict(zip(fields, colors))
    dict_mark = dict(zip(fields, marks))
    dict_ls = dict(zip(['UD', 'DF'], ['solid', 'dashed']))

    for i, vval in enumerate(categ):
        fig, ax = plt.subplots(figsize=(14, 8))
        for field in fields:
            idx = data['target_name'] == field
            sel = data[idx]
            dd_type = sel['DD_type'].unique()[0]
            ax.plot(sel['season'], sel['nvisits_{}'.format(vval)],
                    color=dict_col[field], marker=dict_mark[field],
                    linestyle=dict_ls[dd_type],
                    mfc='None', ms=10, label=field.split(prefix)[-1])

        ax.grid(visible=True)
        ax.set_xlabel(r'season')
        ax.set_ylabel(
            r'$\frac{N_{visits}^{exp}}{N_{visits}^{simu}}$'+value[i]+'1')
        ax.legend(bbox_to_anchor=(1.01, 1.10),
                  ncol=6, fontsize=15, frameon=False)
    plt.show()
