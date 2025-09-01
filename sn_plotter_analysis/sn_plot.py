#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 15:14:35 2025

@author: philippe.gris@clermont.in2p3.fr
"""

import numpy as np
from . import plt
import pandas as pd
from sn_plotter_analysis.sn_analyser_tools import count_all
from sn_plotter_analysis.sn_analyser_tools import clean_level


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


def plot_nsn_tot(nsn_a, config,
                 cumul=False,
                 yvar='nsn', ylab='$\Sigma N_{SN}$',
                 yvar_err='', fields=['COSMOS']):
    """
    Function to plot nsn (no sel) vs year

    Parameters
    ----------
    nsn_a : pandas df
        Data to process.
    config : pandas df
        configuration for the plot.

    cols : list(str), optional
        List of cols (groupby) to estimate nsn. The default is ['year', 'dbName'].
    cumul : bool, optional
        To plot cumulative or not. The default is False.
    fields : list(str), optional
        List of DDFs to consider. The default is ['COSMOS'].

    Returns
    -------
    None.

    """

    plot_nsn_year_all(nsn_a, config,
                      xvar='year', xlab='year',
                      yvar=yvar, ylab=ylab, yvar_err=yvar_err,
                      cumul=cumul, figtit=','.join(fields))


def plot_ddf_year(data, config,
                  cols=['year', 'dbName'],
                  fields=['COSMOS', 'CDFS',
                          'XMM-LSS',
                          'ELAISS1', 'EDFS_a', 'EDFS_b']):
    """
    Function to plot nsn vs year

    Parameters
    ----------
    data : pandas df
        Data to process.
    config : pandas df
        configuration for the plot.
    cols : list(str), optional
        columns to select data. The default is ['year', 'dbName']
    fields : list(str), optional
        List of DDFs to consider. The default is
        ['COSMOS', 'CDFS','XMM-LSS','ELAISS1', 'EDFS_a', 'EDFS_b'].

    Returns
    -------
    None.

    """

    idx = data['field'].isin(fields)
    data = pd.DataFrame(data[idx])
    # plot nsn vs year - with stat error

    datab = count_all(
        data, cols, var=['nsn', 'survey_area'], err_var=['err_nsn'])

    datab['nsn_sqdeg'] = datab['nsn']/datab['survey_area']
    datab['err_nsn_sqdeg'] = datab['err_nsn']/datab['survey_area']
    for cumul in [False, True]:
        plot_nsn_tot(datab, config, yvar='nsn', yvar_err='err_nsn',
                     cumul=cumul, fields=fields)

    ylab = '$N_{SN}/deg^2$'
    plot_nsn_tot(datab, config, yvar='nsn_sqdeg',
                 yvar_err='err_nsn_sqdeg', ylab=ylab,
                 cumul=False, fields=fields)

    # zmin > 0.8, sigmac<=0.04
    var = ['nsn_z_08_sigmaC', 'survey_area']
    err_var = 'err_nsn_z_08_sigmaC'
    datab = count_all(data, cols, var=var, err_var=[err_var])
    datab = clean_level(datab)

    ylab_add = '$z \geq $'+'{}'.format(0.8)
    ylab_add += ', $\sigma_C \leq $'+'{}'.format(0.04)
    ylab = '$\Sigma N_{SN}$'
    if ylab_add != '':
        ylab += '({})'.format(ylab_add)

    for cumul in [False, True]:
        plot_nsn_tot(datab, config, yvar=var[0],
                     yvar_err=err_var, ylab=ylab,
                     cumul=cumul, fields=fields)

    datab['nsn_sqdeg'] = datab['nsn_z_08_sigmaC']/datab['survey_area']
    datab['err_nsn_sqdeg'] = datab['err_nsn_z_08_sigmaC']/datab['survey_area']
    ylab = '$N_{SN}/deg^2$'+'({})'.format(ylab_add)
    plot_nsn_tot(datab, config, yvar='nsn_sqdeg',
                 yvar_err='err_nsn_sqdeg', ylab=ylab,
                 cumul=False, fields=fields)

    var = 'nsn_z_08'
    err_var = 'err_nsn_z_08'
    datac = count_all(data, cols, var=[var], err_var=[err_var])
    datac = clean_level(datac)
    datab = datab.drop(columns=['survey_area'])

    data_m = datac.merge(datab, left_on=['dbName', 'year'], right_on=[
                         'dbName', 'year'], suffixes=['', ''])

    data_m['ratio'] = data_m['nsn_z_08_sigmaC']/data_m['nsn_z_08']
    vv = data_m['err_nsn_z_08_sigmaC']**2
    vv /= data_m['nsn_z_08']**2
    vvb = (data_m['nsn_z_08_sigmaC']/data_m['nsn_z_08']**2)**2
    vvb *= data_m['err_nsn_z_08']**2
    data_m['err_ratio'] = np.sqrt(vvb+vv)

    ylab = '$\\frac{N_{SN}^{z \geq ' + '{}'.format(0.8)
    ylab += ',\sigma_C \leq '+'{}'.format(0.04)
    ylab += '}}{N_{SN}^{z \geq '+'{}'.format(0.8)+'}}$'

    plot_nsn_tot(data_m, config, yvar='ratio',
                 yvar_err='err_ratio', ylab=ylab,
                 cumul=False, fields=fields)
