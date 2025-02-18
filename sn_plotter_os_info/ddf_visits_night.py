#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 10:01:53 2025

@author: philippe.gris@clermont.in2p3.fr
"""
import numpy as np
from . import plt
import pandas as pd


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

    print(data.columns)

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
