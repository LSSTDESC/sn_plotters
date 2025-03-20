#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 14:29:53 2025

@author: philippe.gris@clermont.in2p3.fr
"""
import numpy as np
import pandas as pd
from . import plt
from sn_tools.sn_io import checkDir


def pixelSize(nside):
    """
    Method to return the pixel size

    Parameters
    ----------
    nside : int
        nside healpix param.

    Returns
    -------
    pixSize : float
        pixel size.

    """

    import healpy as hp
    pixSize = hp.nside2pixarea(nside, degrees=True)

    return pixSize


def plot_versus(df, xvar='year', xlabel='year',
                yvar='nsn', ylabel='$N_{SN}$', fig=None, ax=None,
                label='', ls='solid', marker='o', color='k', cumul=False, mfc='k'):
    """
    Function to plot yvar vs xvar

    Parameters
    ----------
    df : pandas df
        Data to plot.
    xvar : str, optional
        x-axis variable. The default is 'year'.
    xlabel : str, optional
        x-axis label. The default is 'year'.
    yvar : str, optional
        y-axis variable. The default is 'nsn'.
    ylabel : str, optional
        y-axis legend. The default is '$N_{SN}$'.
    fig : matplotlib figure, optional
        plot figure. The default is None.
    ax : matplotlib axis, optional
        plot axis. The default is None.
    label : str, optional
        plot label. The default is ''.
    ls : str, optional
        linestyle. The default is 'solid'.
    marker : str, optional
        marker. The default is 'o'.
    color : str, optional
        color. The default is 'k'.
    cumul : bool, optional
        to plot cumilative. The default is False.
    mfc : str, optional
        marker font color. The default is 'k'.

    Returns
    -------
    None.

    """

    if fig is None:
        fig, ax = plt.subplots(figsize=(12, 8))

    ypl = df[yvar]
    if cumul:
        ypl = np.cumsum(ypl)

    print('plotting here', xvar, yvar)
    ax.plot(df[xvar], ypl, ls=ls, marker=marker,
            color=color, label=label, mfc=mfc, markersize=9, lw=2)


def plotMollview(data, varName, figtit, xmin, xmax,
                 nside=128, outDir='.', saveName=''):
    """
    Function to display results as a Mollweid map

    Parameters
    ---------------
    data: pandas df
      data to consider
    varName: str
      name of the variable to display
    leg: str
      legend of the plot
    xmin: float
      min value for the display
    xmax: float
     max value for the display
    nside: int, optional
        nside parameter for healpix. The default is 128
    saveName:str, optional.
       output name for the jpeg. The default is ''

    """
    import healpy as hp
    npix = hp.nside2npix(nside)

    fig = plt.figure(figsize=(8, 6))

    hpxmap = np.zeros(npix, dtype=float)
    hpxmap = np.full(hpxmap.shape, 0.)
    hpxmap[data['healpixID'].astype(
        int)] += data[varName]

    print(np.where(hpxmap < 0.01))

    norm = plt.cm.colors.Normalize(xmin, xmax)
    cmap = plt.cm.jet
    cmap.set_under('w')

    hp.mollview(hpxmap, fig=fig, min=xmin, max=xmax, cmap=cmap,
                title=figtit, nest=True, norm=norm)
    hp.graticule()

    if saveName != '':
        plt.savefig('{}/{}'.format(outDir, saveName))

    plt.close()


def process_WFD(conf_df, dataType, dbDir_WFD, runType,
                timescale_file, timeslots, norm_factor, fName):
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

    Returns
    -------
    wfd : pandas df
        Output data.

    """

    OS_WFDs = conf_df['dbName_WFD'].unique()
    print('dbNames', OS_WFDs)

    # fig, ax = plt.subplots(figsize=(14, 8))
    from_to_load = 'from sn_plotter_analysis.sn_analyser_tools'
    mod_to_load = '{} import load_{}'.format(from_to_load, dataType)
    exec(mod_to_load)
    for OS_WFD in OS_WFDs:
        idx = conf_df['dbName_WFD'] == OS_WFD
        tt = 'load_{}(\'{}\',\'{}\',\'{}\',\'{}\',{},norm_factor={})'.format(
            dataType, dbDir_WFD, OS_WFD, runType,
            timescale_file, timeslots, norm_factor)
        wfda = eval(tt)
        wfda['dbName'] = OS_WFD
        wfda.to_hdf(fName, key='nsn_WFD')
        del wfda


def plot_summary_wfd(wfda, conf_df, timescale='season',
                     cumul=False):
    """
    Method to plot nsn vs year

    Parameters
    ----------
    wfd: pandas df
        Data to process.
    conf_df: pandas df
        config for plot.
    timescale: str, optional
        Time scale to use(season/year). The default is 'season'.
    cumul: bool, optional
        To plot cumulative results. The default is False.

    Returns
    -------
    None.

    """

    wfd = wfda.groupby(['dbName', timescale])[
        'nsn', 'nsn_cosmo'].sum().reset_index()

    fig, ax = plt.subplots(figsize=(18, 8))
    fig.subplots_adjust(right=0.75)
    for dbName in wfd['dbName'].unique():
        idx = wfd['dbName'] == dbName
        sel = wfd[idx]
        idc = conf_df['dbName_WFD'] == dbName
        selp = conf_df[idc]
        ls = selp['ls'].values[0]
        marker = selp['marker'].values[0]
        color = selp['color'].values[0]
        dbNameb = dbName
        plot_versus(sel, fig=fig, ax=ax, cumul=cumul,
                    ls=ls, marker=marker, color=color, mfc=color, label=dbNameb)
        labelb = dbNameb+' - '+'$\sigma_{\mu}\leq \sigma_{int}$'
        labelb = dbNameb+' - '+'$\sigma_C \leq 0.04$'
        labelb = dbNameb+' - '+'cosmo'
        plot_versus(sel, yvar='nsn_cosmo', fig=fig, ax=ax, cumul=cumul,
                    ls='dotted', marker=marker, color=color,
                    mfc='None', label='')

    ax.grid()
    ax.set_xlim([0.95, 10.05])
    ax.set_xlabel(timescale, fontweight='bold')
    legy = '$N_{SN}$'
    if cumul:
        '$\Sigma N_{SN}$'
    ax.set_ylabel(legy)
    # 0, 1.15 for multiple OS
    # ax.legend(loc='upper left', bbox_to_anchor=(
    #    0.1, 1.1), ncol=3, fontsize=15, frameon=False)
    ax.legend(loc='upper center',
              bbox_to_anchor=(1.20, 0.7),
              ncol=1, fontsize=12, frameon=False)
    if cumul:
        xmin, xmax = ax.get_xlim()

        nsn = 1.e6
        ax.plot([xmin, xmax], [nsn, nsn],
                color='dimgrey', lw=2, linestyle='solid')
        ax.text(5, 1.02e6, '1 million SNe Ia', color='dimgrey', fontsize=12)
        nsn = 200000
        ax.plot([xmin, xmax], [nsn, nsn],
                color='dimgrey', lw=2, linestyle='solid')
        ax.text(5, 0.22e6, '200k SNe Ia', color='dimgrey', fontsize=12)


def plot_mollview_wfd(data, timescale, timeslots, nside, varp='nsn', outDir='.'):
    """
    Function to make Mollweid plots for nsn in the WFD survey

    Parameters
    ----------
    data : pandas df
        Data to process.
    timescale : str
        Time scale (year/season).
    timeslots : list
        List of season/years to plot.
    nside : int
        healpix nside parameter.
    varp : str, optional
        var to plot. The default is 'nsn'.
    outDir : str, optional
        output directory to save the plot. The default is '.'.

    Returns
    -------
    None.

    """

    dbNames = data['dbName'].unique()

    varleg = 'N$_{SN}$='
    for dbName in dbNames:
        idx = data['dbName'] == dbName
        sel = data[idx]
        # plot all seasons
        xmin = sel[varp].min()
        xmax = sel[varp].max()
        nsn = int(np.sum(sel[varp]))
        figtit = '{} \n'.format(dbName)
        figtitb = figtit + varleg
        figtitb += '{}'.format(nsn)
        outDirName = '{}/{}'.format(outDir, dbName)
        checkDir(outDirName)
        saveName = 'nsn.png'
        plotMollview(sel, varp, figtitb, xmin, xmax, nside=nside,
                     outDir=outDirName, saveName=saveName)
        # season by season
        for timesl in timeslots:
            idxb = sel[timescale] == timesl
            selb = sel[idxb]
            xmin = selb[varp].min()
            xmax = selb[varp].max()
            nsn = int(np.sum(selb[varp]))
            figtitb = figtit + '{} {} '.format(timescale, timesl)
            figtitb += varleg+'{}'.format(nsn)
            saveName = 'nsn_{}_{}.png'.format(timescale, timesl)
            plotMollview(selb, varp, figtitb,
                         xmin, xmax, nside=nside,
                         outDir=outDirName, saveName=saveName)


def plot_density_wfd(datam, timescale, timeslots, nside, conf_df,
                     varp='nsn', norm_factor=10, plot_indiv=False):
    """
    Function to plot SN densities

    Parameters
    ----------
    datam : pandas df
        Data to process.
    timescale : str
        Timescale to use.
    timeslots : list(int)
        Time slots to select.
    nside : int
        nside healpix parameter.
    conf_df : pandas df
        Config for the plot.
    varp : str, optional
        Data to consider. The default is 'nsn'.
    norm_factor : float, optional
        WFD norm factor. The default is 10.
    plot_indiv : bool, optional
        to plot indiv. The default is False.

    Returns
    -------
    None.

    """

    print(datam.columns)

    idx = datam[varp] > 0.
    data = datam[idx]
    data[varp] /= norm_factor

    data['healpixID'] = data['healpixID'].astype(int)
    healpixId = data['healpixID'].unique().tolist()
    df_pix = pix_RA_Dec(healpixId, nside)
    data = data.merge(df_pix, left_on=['healpixID'], right_on=[
        'healpixID'], suffixes=['', ''])

    dbNames = data['dbName'].unique()
    ylabel = 'N$_{SN}$/deg$^{2}$'
    if varp == 'nsn_cosmo':
        ylabel = 'N$_{SN}^{cosmo}$/deg$^{2}$'
    fig, ax = plt.subplots(figsize=(14, 8))
    fig.subplots_adjust(right=0.75)
    figb, axb = plt.subplots(figsize=(14, 8))
    figb.subplots_adjust(right=0.75)

    for dbName in dbNames:
        idx = data['dbName'] == dbName
        sel = data[idx]
        idc = conf_df['dbName_WFD'] == dbName
        selp = conf_df[idc]
        ls = selp['ls'].values[0]
        marker = selp['marker'].values[0]
        color = selp['color'].values[0]
        vara = '{}_density_mean'.format(varp)
        varb = '{}_density_std'.format(varp)
        dfa = sel.groupby(['healpixID', 'pixRA', 'pixDec'])[
            varp].sum().reset_index()
        df = get_nsn_dec(dfa, varp, delta_dec=5., nside=nside)
        if plot_indiv:
            plot_density_os_summary(df, vara, varb,
                                    fig=None, ax=None,
                                    ylabel=ylabel, figtit=dbName,
                                    ls=ls, color=color,
                                    marker=marker, label='')

        plot_density_os_summary(df, vara, '', fig=fig, ax=ax,
                                ylabel=ylabel, figtit='',
                                ls=ls, color=color, marker=marker, label=dbName)
        plot_density_os_summary(df, '{}_area'.format(varp), '',
                                fig=figb, ax=axb, ylabel='area [deg$^2$]',
                                figtit='', ls=ls, color=color,
                                marker=marker, label=dbName)

    ax.grid(visible=True)
    ax.set_xlabel(r'Dec [deg]')
    ax.set_ylabel(r'{}'.format(ylabel))
    ax.legend(loc='upper center',
              bbox_to_anchor=(1.2, 0.7),
              ncol=1, fontsize=15, frameon=False)

    axb.grid(visible=True)
    axb.set_xlabel(r'Dec [deg]')
    axb.set_ylabel(r'area [deg$^2$]')
    axb.legend(loc='upper center',
               bbox_to_anchor=(1.2, 0.7),
               ncol=1, fontsize=15, frameon=False)


def plot_density_wfd_season(datam, timescale, timeslots, nside, conf_df,
                            varp='nsn', norm_factor=10):
    """
    Function to plot SN densities

    Parameters
    ----------
    datam : pandas df
        Data to process.
    timescale : str
        Timescale to use.
    timeslots : list(int)
        Time slots to select.
    nside : int
        nside healpix parameter.
    conf_df : pandas df
        Config for the plot.
    varp : str, optional
        Data to consider. The default is 'nsn'.
    norm_factor : float, optional
        WFD norm factor. The default is 10.

    Returns
    -------
    None.

    """

    print(datam.columns)

    idx = datam[varp] > 0.
    data = datam[idx]
    data[varp] /= norm_factor

    data['healpixID'] = data['healpixID'].astype(int)
    healpixId = data['healpixID'].unique().tolist()
    df_pix = pix_RA_Dec(healpixId, nside)
    data = data.merge(df_pix, left_on=['healpixID'], right_on=[
        'healpixID'], suffixes=['', ''])

    dbNames = data['dbName'].unique()
    ylabel = 'N$_{SN}$/deg$^{2}$'
    if varp == 'nsn_cosmo':
        ylabel = 'N$_{SN}^{cosmo}$/deg$^{2}$'

    vara = '{}_density_mean'.format(varp)
    varb = '{}_density_std'.format(varp)
    years = range(1, 13)
    marks = ['o', 's', 'P', 'v', '^', 'p']*2
    lsb = ['solid']*6+['dashed']*6
    clrs = ['k']*2+['r']*2+['g']*2+['b']*2+['m']*2+['orange']*2

    mm = dict(zip(years, marks))
    ll = dict(zip(years, lsb))
    ccol = dict(zip(years, clrs))

    for dbName in dbNames:
        idx = data['dbName'] == dbName
        sel = data[idx]
        idc = conf_df['dbName_WFD'] == dbName
        selp = conf_df[idc]

        fig, ax = plt.subplots(figsize=(12, 8))
        fig.subplots_adjust(right=0.85)

        sel = sel.sort_values(by=[timescale])
        slots = sel[timescale].unique()
        for slot in slots:
            idxb = sel[timescale] == slot
            selb = sel[idxb]
            dfa = selb.groupby(['healpixID', 'pixRA', 'pixDec'])[
                varp].sum().reset_index()
            df = get_nsn_dec(dfa, varp, delta_dec=5., nside=nside)
            ls = ll[slot]
            color = ccol[slot]
            marker = mm[slot]

            plot_density_os_summary(df, vara, varb,
                                    fig=fig, ax=ax,
                                    ylabel=ylabel, figtit=dbName,
                                    ls=ls, color=color,
                                    marker=marker,
                                    label='{} {}'.format(timescale, slot))

        ax.grid(visible=True)
        ax.set_xlabel(r'Dec [deg]')
        ax.set_ylabel(r'{}'.format(ylabel))
        ax.legend(loc='upper center',
                  bbox_to_anchor=(1.12, 0.7),
                  ncol=1, fontsize=15, frameon=False)


def plot_density_os_summary(df, varm, varstd, fig=None, ax=None,
                            ylabel='N$_{SN}$/deg$^{2}$',
                            figtit='', ls='None', color='k', marker='o',
                            label=''):
    """
    Function to make the plot

    Parameters
    ----------
    df : pandas df
        Data to plot.
    varm : str
        y-axis var mean.
    varstd : str
        y-axis vzr std.
    fig : matplotlib figure, optional
        Figure for the plot. The default is None.
    ax : matplotlib axis, optional
        axis for the plot. The default is None.
    ylabel : str, optional
        y-label. The default is 'N$_{SN}$/deg$^{2}$'.
    figtit : str, optional
        Figure title. The default is ''.
    ls : str, optional
        Line style. The default is 'None'.
    color : color, optional
        color for the plot. The default is 'k'.
    marker : str, optional
        marker for the plot. The default is 'o'.
    label : str, optional
        label for legend. The default is ''.

    Returns
    -------
    None.

    """

    draw_indiv = False
    if fig is None:
        draw_indiv = True
        fig, ax = plt.subplots(figsize=(12, 8))
    fig.suptitle(figtit)

    ax.plot(df['dec'], df['{}'.format(varm)], color=color,
            linestyle=ls, marker=marker, mfc='None', label=label, markersize=10)
    # ax.errorbar(df['dec'], df['nsn_density_mean'],
    #            yerr=df['nsn_density_std'], color='k')
    if draw_indiv:
        df['plus'] = df['{}'.format(varm)]+df['{}'.format(varstd)]
        df['minus'] = df['{}'.format(varm)]-df['{}'.format(varstd)]
        ax.fill_between(df['dec'], df['plus'],
                        df['minus'], color='yellow')

        ax.grid(visible=True)
        ax.set_xlabel(r'Dec [deg]')
        ax.set_ylabel(r'{}'.format(ylabel))


def pix_RA_Dec(healpixId, nside):
    """
    Function to grab (pixRA,pixDec) from healpixID list

    Parameters
    ----------
    healpixId : list(int)
        List of healpixIDs.
    nside : int
        nside healpix parameter.

    Returns
    -------
    data : pandas df
        results: df['healpixID','pixRA','pixDec'].

    """

    import healpy as hp
    coord = hp.pix2ang(nside, healpixId, nest=True, lonlat=True)
    df_pix = pd.DataFrame(healpixId, columns=['healpixID'])

    df_pix['pixRA'] = coord[0]
    df_pix['pixDec'] = coord[1]

    return df_pix


def get_nsn_dec(data, varp='nsn', delta_dec=5., nside=64):
    """
    Function to estimate the varp density and area per Dec slices

    Parameters
    ----------
    data : pandas df
        Data to process.
    varp : str, optional
        var to consider. The default is 'nsn'.
    delta_dec : float, optional
        dec slice width. The default is 5..
    nside: int, optional
        nside healpix parameter. The default is 64.

    Returns
    -------
    df : pandas df
        output data.

    """

    decs = np.arange(-80., 20., delta_dec)
    bin_centers = (decs[: -1] + decs[1:])/2
    df = pd.DataFrame(bin_centers, columns=['dec'])
    df['dec'] -= delta_dec/2.

    group = data.groupby(pd.cut(data['pixDec'], decs))

    pixSize = pixelSize(nside)
    df[f'{varp}_sum'] = group[varp].sum().to_list()
    df[f'{varp}_density_mean'] = group[varp].mean().to_list()
    df[f'{varp}_density_std'] = group[varp].std().to_list()
    df[f'{varp}_density_mean'] /= pixSize
    df[f'{varp}_density_std'] /= pixSize
    df[f'{varp}_area'] = group.size().to_list()
    df[f'{varp}_area'] *= pixSize

    return df
