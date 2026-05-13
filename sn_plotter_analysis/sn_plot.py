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
import warnings
warnings.filterwarnings("ignore")


def plot_nsn_year_all(nsn, config,
                      xvar='year', xlab='year',
                      yvar='nsn', ylab='N$_{SN}$',
                      yvar_err='',
                      cumul=False, figtit='', os_ref='None'):
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
    os_ref: str, optional
       Ref os to normalize the results. The default is 'None'

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots(figsize=(16, 8))
    fig.subplots_adjust(right=0.78)

    if os_ref != 'None':
        figtit += '\n ref:{}'.format(os_ref)
    fig.suptitle(figtit)

    ref_os = pd.DataFrame()
    if os_ref != 'None':
        idx = nsn['dbName'] == os_ref
        ref_os = nsn[idx]
        if cumul:
            ref_os = get_cumul(ref_os, xvar, yvar, yvar_err)

    config = config.sort_values(by=['dbName_plot'])
    dbNames = config['dbName_plot'].unique()

    for dbName in dbNames:
        if os_ref == dbName:
            continue
        idx = nsn['dbName'] == dbName
        sel = nsn[idx]

        if cumul:
            sel = get_cumul(sel, xvar, yvar, yvar_err)

        # print(sel[['dbName', xvar, yvar, yvar_err]])
        sel = get_norm(sel, ref_os, xvar, yvar, yvar_err)

        toplot = sel[yvar]
        yerr = None
        if yvar_err != '':
            yerr = sel[yvar_err]
        """
        if cumul:
            toplot = np.cumsum(toplot)
            if yerr is not None:
                yerr = np.sqrt(np.cumsum(yerr**2))
        """
        # get config for plot
        idxb = config['dbName_plot'] == dbName
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


def get_sum(data, os_ref, timescale, yvar):
    """
    Function to estimate the correction to be applied from a set of OS
    (typical to account for the weather)

    Parameters
    ----------
    data : pandas df
        Data to process.
    os_ref : str
        ref OS.
    timescale : str
        timescale (year/season).
    yvar : str
        variable to estimate.

    Returns
    -------
    pandas df
       output results
    """

    idx = data['dbName'] == os_ref
    df_ref = pd.DataFrame(data[idx])

    tt = data.groupby([timescale]).apply(
        lambda x: get_vals(x, yvar), include_groups=False).reset_index()

    tt = tt.merge(df_ref, left_on=[timescale],
                  right_on=[timescale])

    return tt


def get_vals(grp, yvar):
    """
    Function to get (mean,rms)

    Parameters
    ----------
    grp : pandas df
        Data to process.
    yvar : str
        col of interest.

    Returns
    -------
    res : TYPE
        DESCRIPTION.

    """

    ll = grp[yvar].to_list()
    mean = grp[yvar].mean()
    std = grp[yvar].std()

    rr = [(mean, std)]
    cols = ['{}_mean'.format(yvar), '{}_std'.format(yvar)]
    res = pd.DataFrame(rr, columns=cols)

    return res


def get_cumul(dfa, xvar, yvar, yvar_err):
    """
    Function to estimate cumulative

    Parameters
    ----------
    dfa : pandas df
        Data to cumulate.
    xvar : str
        x-axis variable.
    yvar : str
        var to cumulate.
    yvar_err : str
        var_err to cumulate.

    Returns
    -------
    dft : TYPE
        DESCRIPTION.

    """
    if len(dfa) == 0:
        return dfa
    dfa = dfa.sort_values(by=[xvar])

    vva = dfa[xvar].to_list()

    dft = pd.DataFrame(vva, columns=[xvar])
    dft[yvar] = np.cumsum(dfa[yvar]).tolist()
    err = np.array(dfa[yvar_err].to_list())
    dft[yvar_err] = np.sqrt(np.cumsum(err**2))
    dft['dbName'] = dfa['dbName'].to_list()

    return dft


def get_norm(df, ref_os, xvar, yvar, yvar_err):
    """
    Function to normalize the results

    Parameters
    ----------
    df : pandas df
        Data to process.
    ref_os : pandas df
        Data of the ref OS.
    xvar : str
        x-axis var.
    yvar : str
        var to normalize.
    yvar_err : str
        err of the variable to normalize.

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    """

    if len(ref_os) >= 1:
        df = df.merge(ref_os, left_on=[xvar], right_on=[
                      xvar], suffixes=['', '_ref'])

        yvar_ref = '{}_ref'.format(yvar)
        yvar_err_ref = '{}_ref'.format(yvar_err)
        # print('rr', df[[xvar, yvar, yvar_ref, yvar_err, yvar_err_ref]])

        df[yvar] /= df[yvar_ref]
        vva = (df[yvar_err]/df[yvar_ref])**2
        vvb = ((df[yvar_err_ref]*df[yvar])/(df[yvar_ref])**2)**2
        df[yvar_err] = np.sqrt(vva+vvb)
        # print('bb', df[[xvar, yvar, yvar_ref, yvar_err, yvar_err_ref]])

    return df


def plot_nsn_tot(nsn_a, config,
                 cumul=False,
                 yvar='nsn', ylab='$\Sigma N_{SN}$',
                 yvar_err='', fields=['COSMOS'], os_ref='None'):
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
    os_ref: str, optional.
        ref OS to normalize the results. The default is 'None'.

    Returns
    -------
    None.

    """

    plot_nsn_year_all(nsn_a, config,
                      xvar='year', xlab='year',
                      yvar=yvar, ylab=ylab, yvar_err=yvar_err,
                      cumul=cumul, figtit=','.join(fields), os_ref=os_ref)


def plot_ddf_area(data, config,
                  cols=['year', 'dbName', 'field'],
                  fields=['COSMOS', 'CDFS',
                          'XMM-LSS',
                          'ELAISS1', 'EDFS_a', 'EDFS_b']):
    idx = data['field'].isin(fields)
    data = pd.DataFrame(data[idx])
    # plot nsn vs year - with stat error

    datab = count_all(
        data, cols, var=['nsn', 'survey_area'], err_var=['err_nsn'])

    plot_nsn_tot(datab, config, yvar='survey_area',
                 ylab='survey area [deg2]',
                 yvar_err='',
                 cumul=False, fields=fields)


def plot_ddf_year(data, config,
                  cols=['year', 'dbName'],
                  fields=['COSMOS', 'CDFS',
                          'XMM-LSS',
                          'ELAISS1', 'EDFS_a', 'EDFS_b'], os_ref='None'):
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
    os_ref: str, optional.
        ref OS to normalize the results. The default is 'None'.
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
    ylab = '$\Sigma N_{SN}$'
    if os_ref != 'None':
        ylab = '$\\frac{\Sigma N_{SN}}{\Sigma N_{SN}^{ref}}$'

    for cumul in [False, True]:
        plot_nsn_tot(datab, config, yvar='nsn', yvar_err='err_nsn',
                     ylab=ylab, cumul=cumul, fields=fields, os_ref=os_ref)

    """
    ylab = '$N_{SN}/deg^2$'
    plot_nsn_tot(datab, config, yvar='nsn_sqdeg',
                 yvar_err='err_nsn_sqdeg', ylab=ylab,
                 cumul=False, fields=fields)
    """
    # zmin > 0.8, sigmac<=0.04
    var = ['nsn_z_08_sigmaC', 'survey_area']
    err_var = 'err_nsn_z_08_sigmaC'
    datab = count_all(data, cols, var=var, err_var=[err_var])
    datab = clean_level(datab)

    ylab_add = '$z \geq $'+'{}'.format(0.8)
    ylab_add += ', $\sigma_C \leq $'+'{}'.format(0.04)
    # ylab = '$\Sigma N_{SN}$'
    if ylab_add != '':
        ylab += '({})'.format(ylab_add)

    for cumul in [False, True]:
        plot_nsn_tot(datab, config, yvar=var[0],
                     yvar_err=err_var, ylab=ylab,
                     cumul=cumul, fields=fields, os_ref=os_ref)

    """
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
    """


def get_weather_impact(data, os_ref, xvar='year', yvar='nsn', fields=['COSMOS']):
    """
    Function to estimate the impact of the weather on nsn,err_nsn

    Parameters
    ----------
    data : pandas df
        Data to process.
    os_ref : str
        ref OS.
    xvar : str, optional
        x-axis var. The default is 'year'.
    yvar : str, optional
        y-axis var. The default is 'nsn'.
    fields : list(str), optional
        List of fields to consider. The default is ['COSMOS'].

    Returns
    -------
    None.

    """

    idx = data['field'].isin(fields)
    data = data[idx]

    cols = ['year', 'dbName']

    datab = count_all(data, cols, var=['nsn'], err_var=['err_nsn'])

    tt = get_sum(datab, os_ref, xvar, yvar)

    tt['rat'] = tt['nsn_mean']/tt['nsn']
    tt['rat1'] = tt['nsn_std']/tt['nsn_mean']
    tt['rat2'] = tt['err_nsn']/tt['nsn']
    print('Fields', fields)
    print(tt[['nsn_mean', 'nsn', 'nsn_std', 'err_nsn', 'rat', 'rat1', 'rat2']])
