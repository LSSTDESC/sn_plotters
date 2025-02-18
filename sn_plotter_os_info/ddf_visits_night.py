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

    print(data)

    res_stat = data.groupby(['target_name', 'season']).apply(
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


def plot_stat(data):

    fields = data['target_name'].unique()
    for vval in ['perfect', 'missing', 'excess']:
        fig, ax = plt.subplots(figsize=(14, 8))
        for field in fields:
            idx = data['target_name'] == field
            sel = data[idx]
            ax.plot(sel['season'], sel['nvisits_{}'.format(vval)])

        ax.grid(visible=True)

    plt.show()
