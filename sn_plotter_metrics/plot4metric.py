#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 15:14:11 2023

@author: philippe.gris@clermont.in2p3.fr
"""
from sn_plotter_metrics import plt
import numpy as np

def plot_vs_OS(data, varx='family', vary='time_budget', legy='Time Budget [%]', title='', fig=None, ax=None, label='', color='k', marker='.', ls='solid', mfc='k'):
    """
    Function to plot results vs OS name

    Parameters
    ----------
    data : pandas df
        data to process
    varx : str, optional
        x-axis col value. The default is 'family'.
    vary : str, optional
        y-axis col value. The default is 'time_budget'.
    legy : str, optional
        y-axis legend. The default is 'Time Budget [%]'.
    title : str, optional
        figure title. The default is ''.
    fig : matplotlib figure, optional
        figure where to plot. The default is None.
    ax : matplotlib axis, optional
        axis where to plot. The default is None.
    label : str, optional
        plot label. The default is ''.
    color : str, optional
        color for the plot. The default is 'k'.
    marker : str, optional
        marker for the plot. The default is '.'.
    ls : str, optional
        linestyle for the plot. The default is 'solid'.
    mfc : str, optional
        marker font color for the plot. The default is 'k'.

    Returns
    -------
    None.

    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))

    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.20)

    # ll = ''
    # if label != '':
    #    ll = data['field'].unique()
    ax.plot(data[varx], data[vary], color=color,
            marker=marker, label='{}'.format(label), linestyle=ls, mfc=mfc)

    ax.grid()
    ax.tick_params(axis='x', labelrotation=20., labelsize=12)
    for tick in ax.xaxis.get_majorticklabels():
        tick.set_horizontalalignment("right")

    if label != '':
        ax.legend()
    ax.set_ylabel(legy)


def plot_vs_OS_dual(data, varx='family', vary=['time_budget'], legy=['Time Budget [%]'], title='', fig=None, ax=None, color='k', marker='.', ls='solid'):
    """
    Function to plot two results vs OS name

    Parameters
    ----------
    data : pandas df
        data to process
    varx : str, optional
        x-axis col value. The default is 'family'.
    vary : str, optional
        y-axis col value. The default is 'time_budget'.
    legy : str, optional
        y-axis legend. The default is 'Time Budget [%]'.
    title : str, optional
        figure title. The default is ''.
    fig : matplotlib figure, optional
        figure where to plot. The default is None.
    ax : matplotlib axis, optional
        axis where to plot. The default is None.
    label : str, optional
        plot label. The default is ''.
    color : str, optional
        color for the plot. The default is 'k'.
    marker : str, optional
        marker for the plot. The default is '.'.
    ls : str, optional
        linestyle for the plot. The default is 'solid'.
    mfc : str, optional
        marker font color for the plot. The default is 'k'.

    Returns
    -------
    None.

    """
   
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8), ncols=1, nrows=len(vary))

    fig.suptitle(title)
    fig.subplots_adjust(bottom=0.150, hspace=0.02)

    lsize = 17
    for io, vv in enumerate(vary):
        ax[io].plot(data[varx], data[vv], color=color,
                    marker=marker)
        ax[io].grid()
        ax[io].tick_params(axis='x', labelrotation=20., labelsize=lsize)
        for tick in ax[io].xaxis.get_majorticklabels():
            tick.set_horizontalalignment("right")

        ax[io].set_ylabel(legy[io], size=lsize)
        ax[io].tick_params(axis='y', labelsize=lsize)
        if io == 0:
            ax[io].get_xaxis().set_ticklabels([])

    #plotDir = '../../Bureau/ddf_fbs_2.1'
    #plotName = '{}/cad_sl_{}.png'.format(plotDir, title)
    #fig.savefig(plotName)
    plt.show()


def plot_hist_OS(data, by='family', what='cadence'):
    """
    Function to plot a set of histograms

    Parameters
    ----------
    data : pandas df
        Data to process.
    by : str, optional
        index to plot. The default is 'family'.
    what : str, optional
        Column to plot. The default is 'cadence'.

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots()
    for fam in np.unique(data['family']):
        idx = data['family'] == fam
        sel = data[idx]
        sel = sel.dropna()
        ax.hist(sel[what], histtype='step', bins=range(1, 25))
        ll = sel[what].to_list()
        p = np.percentile(ll, 50)
        pm = np.percentile(ll, 16)
        pp = np.percentile(ll, 84)
        print(fam, len(np.unique(sel['dbName'])), pm, p, pp)


def plot_series(df, title='', varx='family', what=['time_budget', 'field'], leg=['Time budget [%]', 'DD Field']):
    """
    Function to plot a serie of figures

    Parameters
    ----------
    df : pandas df
        data to plot
    title : str, optional
        Figure title. The default is ''.
    varx : str, optional
        x-axis column. The default is 'family'.
    what : str, optional
        y-axis column. The default is ['time_budget', 'field'].
    leg : str, optional
        y-axis legend. The default is ['Time budget [%]', 'DD Field'].

    Returns
    -------
    None.

    """

    for i, vv in enumerate(what):
        plot_vs_OS(df, varx=varx, vary=vv, legy=leg[i], title=title)


def plot_series_fields(df, title='', varx='family', what=['time_budget_field', 'time_budget_rel'], leg=['Field Time budget [%]', 'Relative Field Time budget [%]']):
    """
    Function to plot a serie of plots per field

    Parameters
    ----------
    df : pandas df
        data to plot.
    title : str, optional
        Figure title. The default is ''.
    varx : str, optional
        x-axis column. The default is 'family'.
    what : list(str), optional
        List of y-axis columns to plot. The default is ['time_budget_field', 'time_budget_rel'].
    leg : list(str), optional
        List of y-axis legends. The default is ['Field Time budget [%]', 'Relative Field Time budget [%]'].

    Returns
    -------
    None.

    """

    ls = ['solid', 'dotted', 'dashed', 'dashdot']*2
    marker = ['.', 's', 'o', '^', 'P', 'h']
    colors = ['k', 'r', 'b', 'm', 'g', 'c']
    for i, vv in enumerate(what):
        fig, ax = plt.subplots(figsize=(12, 8))
        fig.subplots_adjust(top=0.90)
        for io, field in enumerate(np.unique(df['field'])):
            idx = df['field'] == field
            sel = df[idx]
            plot_vs_OS(sel, varx=varx, vary=vv,
                       legy=leg[i], title=title, fig=fig, ax=ax, ls=ls[io], label='{}'.format(field), marker=marker[io], mfc='None', color=colors[io])
        ax.legend(bbox_to_anchor=(0.5, 1.17), ncol=3,
                  frameon=False, loc='upper center')
        ax.grid()
        
    plt.show()


def plot_series_median(df, title='', varx='family', what=['time_budget', 'field'], leg=['Time budget [%]', 'DD Field']):
    """
    Function to plot a set of figures with median values

    Parameters
    ----------
    df : pandas df
        data to plot
    title : str, optional
        Figure title. The default is ''.
    varx : str, optional
        x-axis column. The default is 'family'.
    what : str, optional
        y-axis column. The default is ['time_budget', 'field'].
    leg : str, optional
        y-axis legend. The default is ['Time budget [%]', 'DD Field'].

    Returns
    -------
    None.

    """

    df = df.groupby(varx)[what].median().reset_index()
    for i, vv in enumerate(what):
        plot_vs_OS(df, varx=varx, vary=vv, legy=leg[i], title=title)


def plot_series_median_fields(df, title='', varx='family', what=['time_budget', 'field'], leg=['Time budget [%]', 'DD Field']):
    """
    
    Function to plot a serie of plots per field - median values

    Parameters
    ----------
    df : pandas df
        data to plot.
    title : str, optional
        Figure title. The default is ''.
    varx : str, optional
        x-axis column. The default is 'family'.
    what : list(str), optional
        List of y-axis columns to plot. The default is ['time_budget_field', 'time_budget_rel'].
    leg : list(str), optional
        List of y-axis legends. The default is ['Field Time budget [%]', 'Relative Field Time budget [%]'].
   

    Returns
    -------
    None.

    """

    df = df.groupby([varx, 'field'])[what].median().reset_index()

    for io, field in enumerate(np.unique(df['field'])):
        idx = df['field'] == field
        sel = df[idx]
        plot_vs_OS_dual(sel, varx=varx, vary=what,
                        legy=leg, title=field)
    # ax.grid()


def plot_night(df, dbName, field):
    """
    Function to plot fraction of filter allocation

    Parameters
    ----------
    df : pandas df
        Data to plot.
    dbName : str
        db to plot
    field : str
        field name

    Returns
    -------
    None.

    """
    print(df.columns)
    idx = df['dbName'] == dbName
    idx &= df['field'] == field
    sel = df[idx]

    for season in sel['season'].unique():
        ids = sel['season'] == season
        sels = sel[ids]
        fig, ax = plt.subplots(figsize=(12,8))
        fig.suptitle('{} - Season {}'.format(field,season))
        r = sels['filter_alloc'].tolist()[0]
        rb = sels['filter_frac'].tolist()[0]
        rt = []
        for i, val in enumerate(r):
            rt.append((val, rb[i]))
        tab = np.rec.fromrecords(rt, names=['filter_alloc', 'filter_frac'])
        print('hhh', np.sum(tab['filter_frac']))
        idx = tab['filter_frac'] > 0.02
        tab = tab[idx]
        ax.plot(tab['filter_alloc'], tab['filter_frac'])
        ax.set_ylabel('Fraction')
        ax.grid()
        ax.tick_params(axis='x', labelrotation=20., labelsize=12)
        for tick in ax.xaxis.get_majorticklabels():
            tick.set_horizontalalignment("right")
    plt.show()


def plot_indiv(data, dbName, fig=None, ax=None, xvars=['season', 'season'], xlab=['Season', 'Season'], yvars=['season_length', 'cadence_mean'], ylab=['Season length [days]', 'Mean Cadence [days]'], label='', color='k', marker='.', mfc='k'):
    """
    function to plot parameters corresponding to a field.

    Parameters
    ----------
    data : pandas df
        Data to plot
    dbName : str
        db to plot
    fig : matplotlib figure, optional
        Figure of the plot. The default is None.
    ax : matplotlib axis, optional
        axis of the plot. The default is None.
    xvars : list(str), optional
        x-axis columns to plot. The default is ['season', 'season'].
    xlab : list(str), optional
        x-axis legends. The default is ['Season', 'Season'].
    yvars : list(str), optional
        y-axis columns to plot. The default is ['season_length', 'cadence_mean'].
    ylab : list(str), optional
        y-axis legends. The default is ['Season length [days]', 'Mean Cadence [days]'].
    label : str, optional
        label. The default is ''.
    color : str, optional
        color of the plot. The default is 'k'.
    marker : str, optional
        marker. The default is '.'.
    mfc : str, optional
        marker font color. The default is 'k'.

    Returns
    -------
    None.

    """
    

    if fig is None:
        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8))
        fig.suptitle('{} pointings'.format(field))
        fig.subplots_adjust(hspace=0.02)

    # dbName = sel['dbName'].unique()[0]
    # field = sel['field'].unique()[0]
    # fig.suptitle('{} \n {} pointings'.format(dbName, field))

    for io, vv in enumerate(xvars):
        ax[io].plot(data[vv], data[yvars[io]], label=label,
                    marker=marker, mfc=mfc, color=color)
        ax[io].set_ylabel(ylab[io])
        if io == 0:
            ax[io].get_xaxis().set_ticklabels([])
        if io == 1:
            ax[1].set_xlabel('Season')
        ax[io].grid()


def plot_field(df, xvars=['season', 'season'], xlab=['Season', 'Season'], yvars=['season_length', 'cadence_mean'], ylab=['Season length [days]', 'Mean Cadence [days]'], title=''):
    """
    function to plot parameters corresponding to a field  - one plot per db  

    Parameters
    ----------
    df : pandas df
        Data to plot.
    xvars : list(str), optional
        List of x-axis params to plot. The default is ['season', 'season'].
    xlab : list(str), optional
        List of x-axis label. The default is ['Season', 'Season'].
    yvars : list(str), optional
        List of y-axis parameters to plot. The default is ['season_length', 'cadence_mean'].
    ylab : str, optional
        List of y-axis label. The default is ['Season length [days]', 'Mean Cadence [days]'].
    title : str, optional
        Figure title. The default is ''.

    Returns
    -------
    None.

    """
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(14, 8))
    fig.suptitle(title)
    fig.subplots_adjust(hspace=0.02, right=0.75)
    
    for dbName in df['dbName'].unique():
        idx = df['dbName'] == dbName
        sel = df[idx]
        family = sel['family'].unique()[0]
        marker = sel['marker'].unique()[0]
        color = sel['color'].unique()[0]
        plot_indiv(sel, dbName, fig=fig, ax=ax, xvars=xvars, xlab=xlab, yvars=yvars, ylab=ylab,
                   label=family, marker=marker, color=color, mfc='None')

    ax[1].legend(bbox_to_anchor=(1., 1.), ncol=1, frameon=False)

    ax[0].grid()
    ax[1].grid()
    plt.show()

def plot_filter_alloc(flat, family, field):
    """
    Function to plot filter allocation

    Parameters
    ----------
    flat : pandas df
        Data to plot.
    family : str
        Name of the family to consider.
    field : atr
        Name of the field to plot.

    Returns
    -------
    None.

    """
    toplot = ['filter_frac']
    leg = ['Median obs. night frac']
    idx = flat['family'] == family
    idx &= flat['field'] == field
    idx &= flat['filter_frac'] > 0.05
    # idx &= np.abs(flat['season']-1) < 1.e-5
    sel = flat[idx]

    tit = '{} - {}'.format(family, field)
    plot_series(sel, title=tit, varx='filter_alloc', what=toplot, leg=leg)


def plot_cumsum(selb, title='', xvar='zcomp', xleg='$z_{complete}$',
                yvar='nsn', yleg='$N_{SN}$', ascending=False):
    """
    Function to plot the cumulative sum of yvar as a funtion of x-var

    Parameters
    ----------
    selb : pandas df
        data to process
    title : str, optional
        Figure title. The default is ''.
    xvar : str, optional
        x-axis var. The default is 'zcomp'.
    xleg : str, optional
        x-axis legend. The default is '$z_{complete}$'.
    yvar : str, optional
        y-axis variable. The default is 'nsn'.
    yleg : str, optional
        y-axis legend. The default is '$N_{SN}$'.
    ascending : bool, optional
        x-axis results displayed in ascending mode. The default is False.

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots(figsize=(14, 8))
    fig.suptitle(title)
    fig.subplots_adjust(right=0.75)
    from scipy.interpolate import interp1d
    # fig.suptitle(title)
    fig.subplots_adjust(bottom=0.20)
    for dbName in selb['dbName'].unique():
        idx = selb['dbName'] == dbName
        selp = selb[idx]
        family = selp['family'].unique()[0]
        marker = selp['marker'].unique()[0]
        color = selp['color'].unique()[0]
        selp = selp.sort_values(by=[xvar], ascending=ascending)
        cumulnorm = np.cumsum(selp[yvar])/np.sum(selp[yvar])
        ax.plot(selp[xvar], cumulnorm, marker=marker,
                color=color, mfc='None', label=family)
        interp = interp1d(
            cumulnorm, selp[xvar], bounds_error=False, fill_value=0.)
        zcomp = interp(0.95)
        io = selp[xvar] >= zcomp
        print('zcomp', np.median(selp[io][xvar]))
    ax.grid()
    ax.invert_xaxis()
    ax.legend(bbox_to_anchor=(1.4, 0.8), ncol=1, frameon=False)
    ax.set_xlabel(xleg)
    ax.set_ylabel(yleg)