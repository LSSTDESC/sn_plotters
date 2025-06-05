#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 15:14:35 2025

@author: philippe.gris@clermont.in2p3.fr
"""

import numpy as np
from . import plt


def plot_nsn_year_all(nsn, config,
                      xvar='year', xlab='year',
                      yvar='nsn', ylab='N$_{SN}$',
                      yvar_err='',
                      cumul=False, figtit=''):
    """
    main plot

    Parameters
    ----------
    nsn : pandas df
        Data to plot.
    config : pandas df
        config for the plot.
    xvar : str, optional
        x-axis variable. The default is 'year'.
    xlab : str, optional
        x-axis label. The default is 'year'.
    yvar : str, optional
        y-axis variable. The default is 'nsn'.
    ylab : str, optional
        y-axis label. The default is 'N$_{SN}$'.
    cumul : bool, optional
        To plot cumulative. The default is False.
    figtit : str, optional
        Figure title. The default is ''.

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots(figsize=(15, 8))
    fig.subplots_adjust(right=0.78)
    fig.suptitle(figtit)

    dbNames = nsn['dbName'].unique()
    for dbName in dbNames:
        idx = nsn['dbName'] == dbName
        sel = nsn[idx]
        toplot = sel[yvar]
        yerr = None
        if yvar_err != '':
            yerr = sel[yvar_err]
        if cumul:
            toplot = np.cumsum(toplot)
            if yerr is not None:
                yerr = np.sqrt(np.cumsum(yerr**2))
        # get config for plot
        idxb = config['dbName'] == dbName
        selconf = config[idxb]
        ls = selconf['ls'].values[0]
        color = selconf['color'].values[0]
        mark = selconf['marker'].values[0]
        name = selconf['dbName_plot'].values[0]
        """
        ax.plot(sel[xvar], toplot, color=color,
                marker=mark, linestyle=ls, label=name, mfc='None', lw=2, ms=10)
        """
        ax.errorbar(sel[xvar], toplot, yerr=yerr, color=color,
                    marker=mark, linestyle=ls, label=name, mfc='None', lw=2, ms=10)

    ax.grid(visible=True)
    ax.set_xlabel(r'{}'.format(xlab))
    ax.set_ylabel(r'{}'.format(ylab))
    ax.set_xlim([0.9, 10.1])
    ax.legend(loc='center left', bbox_to_anchor=(
        1, 0.5), ncol=1, fontsize=14, frameon=False)
