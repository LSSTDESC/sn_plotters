#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 14:40:29 2023

@author: philippe.gris@clermont.in2p3.fr
"""
import numpy as np
import pandas as pd
from . import plt


def drop_col(df):
    """
    Function to drop level_x columns

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    """
    levels = []
    for i in range(1, 4):
        levels.append('level_{}'.format(i))

    for vv in levels:
        if vv in df.columns:
            df = df.drop([vv], axis=1)

    return df


class Process_OS:
    def __init__(self, data, grpCol,
                 udfs=['COSMOS', 'XMM-LSS'],
                 dfs=['CDFS', 'EDFS', 'EDFS', 'ELAISS1']):
        """
        Class to process data.
        The goal is to estimate means and std of variables

        Parameters
        ----------
        data : pandas df
            Data to process.
        grpCol : list(str)
            list of columns for the grouped estimation
            (typically: grpCol=['season','dbName]).

        Returns
        -------
        None.

        """

        self.udfs = udfs
        self.dfs = dfs

        self.res = data.groupby(grpCol).apply(
            lambda x: self.get_mean_std(x)).reset_index()

    def get_mean_std(self, grp):
        """
        Method to estimate means and rms per grp

        Parameters
        ----------
        grp : pandas df
            Data to process.

        Returns
        -------
        grp_a : pandas df
            Output data.

        """

        idx = grp['MoM'] > 0.

        grp_am = grp[idx]
        grp_ba = grp[~idx]

        grp_a = self.calc(grp_am, calcCol='MoM')

        grp_am['sigma_w0'] = 100.*np.sqrt(grp_am['Cov_w0_w0_fit'])
        grp_am['sigma_wa'] = 100.*np.sqrt(grp_am['Cov_wa_wa_fit'])

        grp_a_w0 = self.calc(grp_am, calcCol='sigma_w0')
        grp_a_wa = self.calc(grp_am, calcCol='sigma_wa')

        grp_n = pd.DataFrame()
        """
        if 'MoM_DETF' in grp.columns:
            grp_n = self.calc(grp_am, calcCol='MoM_DETF')

        grp_ba['sigma_w'] = 100.*np.sqrt(grp_ba['Cov_w0_w0_fit'])
        grp_b = self.calc(grp_ba, calcCol='sigma_w')
        """
        dict_field = {}

        field_list = self.udfs+self.dfs+['all_Fields']

        for field in field_list:
            grp_c = self.calc(grp_am, calcCol=field)
            names = ['{}_mean'.format(field), '{}_std'.format(field)]
            grp_c[names] = grp_c[names].astype(int)
            dict_field[field] = grp_c[names]

        cola = grp_a.columns
        """
        colb = grp_b.columns
        set_diff = set(colb) - set(cola)
        diff = list(set_diff)
        grp_a = pd.concat((grp_a, grp_b[diff]), axis=1)
        """
        grp_a = pd.concat((grp_a, grp_n), axis=1)
        grp_a = pd.concat((grp_a, grp_a_w0), axis=1)
        grp_a = pd.concat((grp_a, grp_a_wa), axis=1)

        for key, vals in dict_field.items():
            grp_a = pd.concat((grp_a, vals), axis=1)

        grp_a = drop_col(grp_a)

        return grp_a

    def calc(self, df, calcCol='sigma_w'):
        """
        Method to estimate mean and rms of a column

        Parameters
        ----------
        df : pandas df
            Data to process.
        calcCol : str, optional
            Column of interest. The default is 'sigma_w'.

        Returns
        -------
        res : pandas df
            columns: calCol_mean, calcCol_std.

        """

        var_mean = '{}_mean'.format(calcCol)
        var_std = '{}_std'.format(calcCol)
        res = pd.DataFrame({var_mean: [df[calcCol].mean()],
                            var_std: [df[calcCol].std()]})
        res = res.replace([np.inf, -np.inf, np.nan], 0)

        return res


def cosmo_plot(df,
               varx='season', legx='season',
               vary='MoM', legy='MoM',
               ax=None, ls='solid', marker='.', color='k', leg='', msize=10):
    """
    Function to make a cosmo plot

    Parameters
    ----------
    df : pandas df
        Data to plot.
    varx : str, optional
        x-axis variable. The default is 'season'.
    legx : str, optional
        x-axis legend. The default is 'season'.
    vary : str, optional
        y-axis variable. The default is 'MoM'.
    legy : str, optional
        y-axis legend. The default is 'MoM'.
    prior : str, optional
        prior value. The default is 'noprior'.
    ax : matplotlib axis, optional
        where the plot show be made. The default is None.
    ls : str, optional
        Line style. The default is 'solid'.
    marker : str, optional
        Marker type. The default is '.'.
    color : str, optional
        color. The default is 'k'.
    leg : str, optional
        Label for legend. The default is ''.
    msize : float, optional
        Marker size. The default is 10.

    Returns
    -------
    None.

    """

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))

    vary_m = '{}_mean'.format(vary)
    vary_std = '{}_std'.format(vary)

    ax.errorbar(df[varx], df[vary_m], yerr=df[vary_std],
                ls=ls, marker=marker, color=color,
                label=leg, markersize=msize, mfc='None')

    ax.grid()
    ax.set_xlabel(legx)
    ax.set_ylabel(legy)


def cosmo_four(resdf, timescale='year'):
    """
    Method to make a plot of 4 variables

    Parameters
    ----------
    resdf : pandas df
        data to plot.

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))
    fig.subplots_adjust(hspace=0., wspace=0.)

    varx = timescale
    legx = timescale

    varys = ['sigma_w', 'MoM']
    priors = ['noprior', 'prior']
    leg = dict(zip(varys, [r'$\sigma_w$', r'$MoM$']))

    ipos = dict(zip(varys, [0, 1]))
    jpos = dict(zip(priors, [0, 1]))

    for vary in varys:
        for prior in priors:
            ix = ipos[vary]
            jx = jpos[prior]
            ido = resdf['prior'] == prior
            cosmo_plot(resdf[ido], varx=varx, legx=legx, vary=vary,
                       legy=leg[vary], prior=prior, ax=ax[ix, jx])
            if jx == 1:
                ax[ix, jx].yaxis.tick_right()
                ax[ix, jx].yaxis.set_label_position("right")
            if ix == 0:
                ax[ix, jx].set_title(prior)


def plot_allOS(resdf, config, dataCol='dbName_DD', configCol='dbName',
               prior='prior', varx='year', legx='year', vary='MoM', legy='$MoM$',
               figtitle='with prior', dbNorm=float('0.3'), leg_prefix=''):
    """
    Function to plot all OS on one single plot

    Parameters
    ----------
    resdf : pandas df
        Data to plot.
    config : pandas df
        Configuration (marker, linestyle, ...).
    dataCol: str, optional.
        Data columns for selection. The default is dbName_DD.
    configCol: str, optional.
        Name of the config column. The default is 'dbName'.
    prior : str, optional
        Prior status (prior or noprior). The default is 'prior'.
    varx : str, optional
        x-axis var. The default is 'year'.
    legx : str, optional
        x-axis label. The default is 'year'.
    vary : str, optional
        y-axis var. The default is 'MoM'.
    legy : str, optional
        y-axis label. The default is '$MoM$'.
    figtitle : str, optional
        Figure suptitle. The default is 'with prior.
    dbNorm: str, optional
       db for normalization

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots(figsize=(14, 8))
    fig.subplots_adjust(right=0.75)

    fig.suptitle(figtitle)

    idx = resdf['prior'] == prior
    sela = resdf[idx]

    if dbNorm != '':
        sela = normalize(sela, dbNorm, dataCol, vary, timescale=varx)

    for i, row in config.iterrows():
        idx = sela[dataCol] == row[configCol]
        sel = sela[idx]
        leg = row[configCol]
        if isinstance(leg, str):
            legs = leg.split('_')
            leg = '_'.join(legs[:-1])
        if leg_prefix != '':
            leg = '{}{}'.format(leg_prefix, leg)
        cosmo_plot(sel, varx=varx, legx=legx, vary=vary,
                   legy=legy, ax=ax, ls=row['ls'],
                   marker=row['marker'], color=row['color'], leg=leg)

    ax.grid(visible=True)
    ax.legend(loc='upper center',
              bbox_to_anchor=(1.20, 0.7),
              ncol=1, fontsize=15, frameon=False)
    # ax.grid()


def normalize(sela, dbNorm, dataCol, vary, timescale='year'):

    selb = sela[[timescale, dataCol, '{}_mean'.format(
        vary), '{}_std'.format(vary)]]

    ido = selb[dataCol] == dbNorm
    selnorm = selb[ido]

    selm = selb.merge(selnorm, left_on=[timescale], right_on=[timescale])

    vvm = '{}_mean'.format(vary)
    vvr = '{}_std'.format(vary)
    selm[vvm] = selm['{}_x'.format(vvm)]/selm['{}_y'.format(vvm)]
    selm[vvr] = 0.
    selm[dataCol] = selm['{}_x'.format(dataCol)]

    return selm
