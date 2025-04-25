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
               vary='MoM_mean', legy='MoM',
               vary_std='MoM_std',
               ax=None, ls='solid', marker='.', color='k', leg='',
               msize=10, comment_on_plot='', fill_between=False):
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
    comment_on_plot: str, optional
      to add a comment on the plot. The default is ''
    fill_between: bool, optional
      to fill +-1 sigma area with yeallo. The default is False.

    Returns
    -------
    None.

    """

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))

    yerr = None
    if vary_std != '':
        yerr = df[vary_std]
    ax.errorbar(df[varx], df[vary], yerr=yerr,
                ls=ls, marker=marker, color=color,
                label=leg, markersize=msize, mfc='None', lw=3)

    if fill_between:
        dfb = pd.DataFrame(df)
        vary_plus = '{}_plus_sigma'.format(vary)
        vary_minus = '{}_minus_sigma'.format(vary)

        dfb[vary_plus] = dfb[vary]+dfb[vary_std]
        dfb[vary_minus] = dfb[vary]-dfb[vary_std]
        ax.fill_between(dfb[varx], dfb[vary_plus],
                        dfb[vary_minus], color='yellow')

    ax.grid()
    ax.set_xlabel(legx, fontsize=25)
    ax.set_ylabel(legy, fontsize=25)


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
            sel = resdf[ido]
            cosmo_plot(sel, varx=varx, legx=legx, vary=vary,
                       legy=leg[vary], prior=prior, ax=ax[ix, jx])

            if jx == 1:
                ax[ix, jx].yaxis.tick_right()
                ax[ix, jx].yaxis.set_label_position("right")
            if ix == 0:
                ax[ix, jx].set_title(prior)


def save_data(resdf, prior='prior', varx='year', year_max=6):
    """
    Function to save SMom values in csv+latex

    Parameters
    ----------
    resdf : pandas df
        Data to process.
    prior : str, optional
        prior or not prior. The default is 'prior'.
    varx : str, optional
        Timescale of the data to save. The default is 'year'.
    year_max : int, optional
        timeslot max. The default is 6.

    Returns
    -------
    None.

    """

    idx = resdf['prior'] == prior
    sela = resdf[idx]

    # save SMoM in csv file+latex output
    idxa = sela[varx] == year_max

    selb = sela[idxa]
    # print(selb[['dbName', 'MoM_mean', 'MoM_std']])

    selb = selb.sort_values(by=['MoM_mean'])
    selb[['dbName', 'dbName_plot', 'MoM_mean', 'MoM_std', 'year']].to_csv(
        'smom_final.csv', index=False)

    print_latex(selb)


def plot_allOS(resdf, config, dataCol='dbName_DD', configCol='dbName',
               prior='prior', varx='year', legx='year', vary='MoM', legy='$MoM$',
               vary_std='MoM_std',
               figtitle='with prior', dbNorm=float('0.3'), leg_prefix='',
               comment_on_plot='', fill_between=False, year_max=6):
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
    fill_between: bool, optional
         to fill +-1 sigma area with yeallo. The default is False.
    year_max: int, optional
       last year of the survey. The default is 6.

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots(figsize=(18, 8))
    fig.subplots_adjust(right=0.71)

    fig.suptitle(figtitle, color='b')

    idx = resdf['prior'] == prior
    sela = resdf[idx]

    """
    # save SMoM in csv file+latex output
    idxa = sela[varx] == year_max

    selb = sela[idxa]
    # print(selb[['dbName', 'MoM_mean', 'MoM_std']])

    selb = selb.sort_values(by=['MoM_mean'])
    selb[['dbName', 'MoM_mean', 'MoM_std', 'year']].to_csv(
        'smom_final.csv', index=False)
    
    print_latex(selb)
    """
    if dbNorm != '':
        sela = normalize(sela, dbNorm, dataCol, vary, vary_std, timescale=varx)
        idx = sela['dbName_DD'] == dbNorm
        sela = sela[~idx]

    config = config.sort_values(by=['dbName'])
    for i, row in config.iterrows():
        if row[configCol] == dbNorm:
            continue
        idx = sela[dataCol] == row[configCol]
        sel = sela[idx]
        leg = row[configCol]
        if isinstance(leg, str):
            legs = leg.split('_')
            leg = '_'.join(legs[:-1])
        if leg_prefix != '':
            leg = '{}{}'.format(leg_prefix, leg)
        leg = row['dbName_plot']
        cosmo_plot(sel, varx=varx, legx=legx, vary=vary,
                   legy=legy, vary_std=vary_std, ax=ax, ls=row['ls'],
                   marker=row['marker'], color=row['color'],
                   leg=leg, comment_on_plot=comment_on_plot,
                   fill_between=fill_between)

    ax.grid(visible=True)
    ax.legend(loc='upper center',
              bbox_to_anchor=(1.25, 0.7),
              ncol=1, fontsize=18, frameon=False)
    ymin, ymax = ax.get_ylim()
    ymin = sel[vary].min()
    ax.text(5, 1.2*ymin, comment_on_plot, color='blue', fontsize=18)
    # ax.grid()


def plot_allOS_survey(res_csv='smom_final.csv', dbNorm='baseline_v3.4_10yrs',
                      dataCol='dbName', varx='year', vary='MoM_mean',
                      vary_std='MoM_std'):
    """
    Function to plot SMoM (relative or absolute) for the full survey

    Parameters
    ----------
    res_csv : str, optional
        SMoM values. The default is 'smom_final.csv'.
    dbNorm : str, optional
        Db name for normalization. The default is 'baseline_v3.4_10yrs'.
    dataCol : str, optional
        Db name col. The default is 'dbName'.
    varx : str, optional
        x-axis var. The default is 'year'.
    vary : str, optional
        y-axis var. The default is 'MoM_mean'.
    vary_std : str, optional
        y-axis std var. The default is 'MoM_std'.

    Returns
    -------
    None.

    """

    # get the data
    data = pd.read_csv(res_csv)

    sela = data
    if dbNorm != '':
        sela = normalize(sela, dbNorm, dataCol, vary, vary_std, timescale=varx)
        idx = sela['dbName'] == dbNorm
        sela = sela[~idx]

    fig, ax = plt.subplots(figsize=(18, 8))
    fig.subplots_adjust(bottom=0.20)

    ttit = ''
    if dbNorm != '':
        ttit = 'ref: {} \n'.format(dbNorm)
    ttit += '{} years'.format(int(data[varx].median()))
    fig.suptitle(ttit, color='b')
    sela = sela.sort_values(by=['MoM_mean'], ascending=False)
    """
    pref = '_v3.4_10yrs'
    sela['dbName'] = sela['dbName'].str.split(pref).str[0]
    pref = 'v3.4_10yrs'
    sela['dbName'] = sela['dbName'].str.split(pref).str[0]
    """
    ax.plot(sela['dbName_plot'], sela['MoM_mean'], color='k',
            linestyle='dotted', marker='o', mfc='r', ms=7, lw=3)

    plt.setp(ax.get_xticklabels(), rotation=30,
             ha="right", rotation_mode="anchor", fontsize=12)

    if dbNorm != '':
        legy = '$\\frac{\\Delta SMoM}{SMoM}$ [%]'
    else:
        legy = 'SMoM'
    ax.grid(visible=True)
    ax.set_ylabel(r'{}'.format(legy), fontsize=25)


def normalize(sela, dbNorm, dataCol, vary, vary_std, timescale='year'):
    """


    Parameters
    ----------
    sela : pandas df
        Data to process.
    dbNorm : str
        OS to use as norm.
    dataCol : str
        OS name col.
    vary : str
        var of interest.
    timescale : str, optional
        Time scale (season/year). The default is 'year'.

    Returns
    -------
    selm : pandas df
        Output data.

    """

    vvm = '{}'.format(vary)
    vvr = '{}'.format(vary_std)
    ccols = [timescale, dataCol, vvm, vvr, 'dbName_plot']
    print('ccols', ccols)
    selb = sela[ccols]

    ido = selb[dataCol] == dbNorm
    selnorm = selb[ido]

    selm = selb.merge(selnorm, left_on=[timescale], right_on=[timescale])

    selm[vvm] = 100.*(selm['{}_x'.format(vvm)]-selm['{}_y'.format(vvm)]
                      )/selm['{}_x'.format(vvm)]
    selm[vvr] = 0.
    selm[dataCol] = selm['{}_x'.format(dataCol)]
    selm['dbName_plot'] = selm['dbName_plot_x']

    print('allo', selm.columns)
    return selm


def print_latex(res):
    """
    Function to print result ready for latex

    Parameters
    ----------
    res : pandad df
        Data to process.

    Returns
    -------
    None.

    """

    print('\\begin{table}[!htbp]')
    print('\\begin{center}')
    print(
        '\caption{\smom~values (ascending order) after a 10-season survey. The standard deviation corresponding to \smom~distribution resulting from the fit of 50 random surveys is also quoted.}\label{tab:smom_final}')
    print('\\begin{tabular}{l|c}')
    print('\hline')
    print('\hline')
    print('Observing Strategy & \smom \\\\')
    print('\hline')
    for i, row in res.iterrows():
        dbName = row['dbName']
        dbName = '_'.join(dbName.split('_')[:-1])
        dbName = dbName.replace('_', '\_')
        print('{} & {} $\pm$ {}\\\\'.format(
            dbName, int(np.rint(row['MoM_mean'])), int(np.rint(row['MoM_std']))))
    print('\hline')
    print('\hline')
    print('\end{tabular}')
    print('\end{center}')
    print('\end{table}')
