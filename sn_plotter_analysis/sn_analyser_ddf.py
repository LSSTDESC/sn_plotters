#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 16:50:25 2025

@author: philippe.gris@clermont.in2p3.fr
"""
import numpy as np
import pandas as pd
from sn_analysis import plt
from sn_analysis.sn_calc_plot import bin_it, bin_it_mean, bin_it_effi


def plot_nsn(ax, selb, xvar, yvar, yvar_cut,
             smoothIt, marks, listy, timescale, timeslot, norm_factor):
    """
    Function to plot nsn vs z

    Parameters
    ----------
    ax : matplotlib axis
        plot axis.
    selb : pandas df
        Data to plot.
    xvar : str
        x-axis var.
    yvar : str
        y-axis var.
    yvar_cut : float
        y var selection cut.
    smoothIt : bool
        To smooth (spline) displayed curves.
    marks : dict
        Markers used for the plot.
    listy : dict
        linestyle used for the plot.
    timescale : str
        time scale to use (season/year).
    timeslot : int
        time slot for display.
    norm_factor : int
        Normalization factor.

    Returns
    -------
    None.

    """

    df = bin_it(selb, xvar=xvar, norm_factor=norm_factor,
                bins=np.arange(0.01, 1.12, 0.1))
    # ax.errorbar(df['z'], df['sigma_mu'], yerr=df['sigma_mu_std'])
    if smoothIt:
        from scipy.interpolate import make_interp_spline
        xnew = np.linspace(
            np.min(df[xvar]), np.max(df[xvar]), 100)
        spl = make_interp_spline(
            df[xvar], df[yvar], k=3)  # type: BSpline
        spl_smooth = spl(xnew)

        ax.plot(xnew, spl_smooth, color='k',
                marker=marks[timeslot], ls=listy[timeslot],
                mfc='None', ms=10, markevery=5,
                label='{} {}'.format(timescale, timeslot))

    else:
        ax.plot(df[xvar], df[yvar], color='k',
                marker=marks[timeslot], ls=listy[timeslot],
                mfc='None', ms=10, markevery=5,
                label='{} {}'.format(timescale, timeslot))

    xmin = 0.2
    xmax = 1.08
    ymin = 0.0
    ymax = None
    ax.set_xlim([xmin, xmax])
    # ax.set_ylim([ymin, ymax])


def get_val(var):
    """
    Function to grab values from parser

    Parameters
    ----------
    var : str
        var to process.

    Returns
    -------
    var : list(int)
        Result.

    """
    if '-' in var:
        seas_spl = var.split('-')
        seas_min = int(seas_spl[0])
        seas_max = int(seas_spl[1])
        var = range(seas_min, seas_max+1)
    else:
        var = var.split(',')
        var = list(map(int, var))

    return var


def plot_survey_features(data, field, dbName, norm_factor, config, nside, timescale, timeslots):
    """
    Function to plot survey features related to sn

    Parameters
    ----------
    data : pandas df
        Data to process.
    norm_factor : int
        normalization factor.
    config : dict
        config params.
    nside : int
        nside param (healpix).
    timescale : str
        time scale to use (season/year).

    Returns
    -------
    None.

    """

    # field = 'COSMOS'
    # dbName = 'baseline_v3.4_10yrs'
    # dbName = 'roll_uniform_early_half_mjdp67_v3.4_10yrs'
    # dbName = 'DDF_DESC_0.80_WZ_0.07'

    plot_sn_features(data, field, dbName, timescale, timeslots,
                     yvar='sigma_mu', ylabel='$\\sigma_{\mu}$ [mag]',
                     type_plot='sigma_mu', smoothIt=False)
    plot_sn_features(data, field, dbName, timescale, timeslots,
                     yvar='NSN', ylabel='$N_{SN}$', type_plot='nsn',
                     smoothIt=True, norm_factor=norm_factor)

    plot_sn_features(data, field, dbName, timescale, timeslots, smoothIt=True)

    """
    df = get_zmax_field(data, field, dbName, timescale, zmin=0.7, sigmaC=1.e6)

    fig, ax = plt.subplots()
    ax.plot(df[timescale], df['nsn_zmin'], 'ko')

    plt.show()
    """


def plot_sn_features(data, field, dbName, timescale, timeslots,
                     xvar='z', xlabel='$z$', yvar='sigma_mu',
                     ylabel='$frac^{N_{SN}}_{\sigma_{\mu} \leq \sigma_{int}}$',
                     yvar_cut=0.12, type_plot='effi', smoothIt=False,
                     norm_factor=1):
    """
    Function to plot sn features from survey

    Parameters
    ----------
    data : pandas df
        Data to process.
    field : str
        Field of interest.
    dbName : str
        OS name.
    timescale : str
        Time scale to use (season/year).
    xvar : str, optional
        x-axis var. The default is 'z'.
    xlabel : str, optional
        x-axis label. The default is '$z$'.
    yvar : str, optional
        y-axis var. The default is 'sigma_mu'.
    ylabel : str, optional
        y-axis label.
        The default is '$frac^{N_{SN}}_{\sigma_{\mu} \leq \sigma_{int}}$'.
    yvar_cut : float, optional
        y-axis selection cut. The default is 0.12.
    type_plot : str, optional
        type of plot (sigma_mu, nsn, effi). The default is 'effi'.
    smoothIt : bool, optional
        To smooth (spline) displayed curves. The default is False.
    norm_factor : int, optional
        normalization factor. The default is 1.

    Returns
    -------
    None.

    """

    idx = data['field'] == field
    idx &= data['dbName'] == dbName

    dbNameb = '_'.join(dbName.split('_')[:-1])
    sel = data[idx]

    # for each year: sigmamu vs z
    fig, ax = plt.subplots(figsize=(12, 8))
    fig.subplots_adjust(right=0.82)
    fig.suptitle('{} - {}'.format(dbNameb, field))
    ttimes = range(1, 12)
    lls = ['solid']*4+['dashed']*4+['dotted']*4
    mmarkers = ['o', '*', '^', 'h']*3
    listy = dict(zip(ttimes, lls))
    marks = dict(zip(ttimes, mmarkers))

    for timeslot in timeslots:

        idxb = sel[timescale] == timeslot
        selb = sel[idxb]

        eval('plot_{}(ax, selb, xvar, yvar, yvar_cut, smoothIt,marks, listy, timescale, timeslot, norm_factor)'.format(type_plot))

        """
        if type_plot == 'sigma_mu':
            plot_sigma_mu(ax, selb, xvar, yvar, yvar_cut, smoothIt,
                          marks, listy, timescale, timeslot)
        if type_plot == 'nsn':
            plot_nsn(ax, selb, xvar, yvar, yvar_cut, smoothIt,
                     marks, listy, timescale, timeslot, norm_factor)

        if type_plot == 'effi':
            plot_effi(ax, selb, xvar, yvar, yvar_cut, smoothIt,
                      marks, listy, timescale, timeslot)
        """
    ax.grid(visible=True)
    ax.set_xlabel(r'{}'.format(xlabel))
    ax.set_ylabel(r'{}'.format(ylabel))

    ax.legend(loc='upper center',
              bbox_to_anchor=(1.12, 0.7),
              ncol=1, fontsize=15, frameon=False)


def plot_DDF_nsn(data, norm_factor, config, nside,
                 timescale='year', yleg_add='', cumul=False,
                 plots=['nsn_field_OS', 'nsn_OS', 'pix_area']):
    """


    Parameters
    ----------
    data : pandas df
        Data to plot.
    norm_factor : float
        norm factor.
    config : pandas df
        config for plots.
    nside : int
        nside healpix parameter.
    timescale : str, optional
        time scale for plots. The default is 'year'.
    yleg_add : str, optional
        additionnal y-axis legend. The default is ''.
    cumul : bool, optional
        to display cumulative nsn. The default is False.
    plots : list(str), optional
        List of plots to display. 
        The default is ['nsn_field_OS', 'nsn_OS', 'pixarea'].

    Returns
    -------
    None.

    """

    # mypl = Plot_nsn_vs(data, norm_factor, nside)
    # mypl.plot_nsn_mollview()
    """
    # mypl.plot_nsn_versus_two(xvar='z', xleg='z', logy=True,
    #                         cumul=True, xlim=[0.01, 1.1])
    mypl.plot_nsn_mollview()
    """

    # estimate the number of sn for all the fields/season

    sums = get_sums_nsn(data, norm_factor, nside, cols=[
        timescale, 'dbName', 'field'])
    sumt = get_sums_nsn(data, norm_factor, nside, cols=[timescale, 'dbName'])

    sumb = get_sums_nsn(data, norm_factor, nside, cols=['dbName'])

    print(sumb)
    # plot_field(sums, mypl, config, xvar=timescale,
    #           xleg=timescale, cumul=True)
    # plot_field(sums, mypl, xvar=timescale, xleg=timescale,
    #           yvar='pixArea', yleg='Observed Area [deg$^{2}$]')

    # total number of SN per season/OS
    yleg = '$N_{SN}$'
    if cumul:
        yleg = '$\Sigma N_{SN}$'
    yleg += yleg_add

    if 'nsn_field_OS' in plots:
        plot_field(sums, config, xvar=timescale, xleg=timescale,
                   cumul=cumul, yleg=yleg)

    # total number of SN per season/OS
    if 'nsn_os' in plots:
        plot_field(sumt, config, xvar=timescale, xleg=timescale,
                   cumul=cumul, yleg=yleg)

    if 'pix_area' in plots:
        plot_field(sumt, config, xvar=timescale, xleg=timescale,
                   yvar='pixArea', yleg='Observed Area [deg$^{2}$]')
    # plt.show()


def plot_field(data, config, xvar='season', xleg='season',
               yvar='nsn', yleg='$N_{SN}$', cumul=False, norm='', logy=False):
    """
    Function to plot a set of fields results

    Parameters
    ----------
    data : array
        Data to process.
    config: pandas df
      config for plots
    xvar : str, optional
        x-axis variable. The default is 'season'.
    xleg : str, optional
        x-axis label. The default is 'season'.
    yvar : str, optional
        y-axis var. The default is 'nsn'.
    yleg : str, optional
        y-axis label. The default is '$N_{SN}$'.
    cumul : bool, optional
        for cumulative plot. The default is False.
    Returns
    -------
    None.

    """

    if norm != '':
        # normalize the results here
        idx = data['dbName'] == norm
        selnorm = data[idx]
        vmerge = ['field', xvar]
        df = data.merge(selnorm, left_on=vmerge, right_on=vmerge)
        df['{}'.format(yvar)] = df['{}_x'.format(yvar)]/df['{}_y'.format(yvar)]
        df['dbName'] = df['dbName_x']
        data = pd.DataFrame(df)

    for field in data['field'].unique():
        idx = data['field'] == field
        sela = data[idx]
        fig, ax = plt.subplots(figsize=(14, 8))
        for dbName in sela['dbName'].unique():
            idxb = sela['dbName'] == dbName
            selb = sela[idxb]
            idxc = config['dbName_DD'] == dbName
            conf = config[idxc]
            ls = conf['ls'].to_list()[0]
            color = conf['color'].to_list()[0]
            marker = conf['marker'].to_list()[0]
            plot_versus(selb, xvar, xleg,
                        yvar, yleg,
                        figTitle=field, label=dbName,
                        fig=fig, ax=ax, xlim=None, cumul=cumul,
                        ls=ls, color=color,
                        marker=marker)

        ax.legend()
        # ax.grid()
        ax.set_xlabel(xleg, fontweight='bold')
        ax.set_ylabel(yleg, fontweight='bold')
        ax.grid(visible=True)
        if logy:
            ax.set_yscale("log")


def plot_sigma_mu(ax, selb, xvar, yvar, yvar_cut,
                  smoothIt, marks, listy, timescale, timeslot, norm_factor):
    """
    Function to plot sigma_mu vs z

    Parameters
    ----------
    ax : matplotlib axis
        plot axis.
    selb : pandas df
        Data to plot.
    xvar : str
        x-axis var.
    yvar : str
        y-axis var.
    yvar_cut : float
        y var selection cut.
    smoothIt : bool
        To smooth (spline) displayed curves.
    marks : dict
        Markers used for the plot.
    listy : dict
        linestyle used for the plot.
    timescale : str
        time scale to use (season/year).
    timeslot : int
        time slot for display.
    norm_factor : int
        Normalization factor.
    Returns
    -------
    None.

    """

    df = bin_it_mean(selb, xvar=xvar, yvar=yvar,
                     bins=np.arange(0.01, 1.12, 0.07))
    # ax.errorbar(df['z'], df['sigma_mu'], yerr=df['sigma_mu_std'])
    ax.plot(df[xvar], df[yvar], color='k', marker=marks[timeslot],
            ls=listy[timeslot], mfc='None', ms=10, markevery=5,
            label='{} {}'.format(timescale, timeslot))
    xmin = 0.2
    xmax = 1.08
    ymin = 0.0
    ymax = 0.6
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.plot([xmin, xmax], [0.12]*2, ls='dashed', color='r')
    ttext = '$\sigma_{int}\sim$0.12'
    ax.text(0.3, 0.13, ttext, color='r')


def get_sums_nsn(data, norm_factor, nside, cols=['season', 'dbName', 'field']):
    """
    Function to estimate global parameters (nsn, pixArea, ...)

    Parameters
    ----------
    data : pandas df
        Data to process.
    norm_factor : float
        Normalization factor.
    nside : int
        healpix nside parameter.
    cols : list(str), optional
        Columns to make groups. The default is ['season', 'dbName', 'field'].

    Returns
    -------
    sums : pandas df
        Output data.

    """

    sums = data.groupby(cols).size().to_frame('nsn').reset_index()
    sums['nsn'] /= norm_factor

    if 'field' not in cols:
        sums['field'] = ','.join(data['field'].unique())
    pix = data.groupby(cols).apply(
        lambda x: pd.DataFrame({'npixels': [len(x['healpixID'].unique())]})).reset_index()

    pix['pixArea'] = pixelSize(nside)*pix['npixels']

    sums = sums.merge(pix, left_on=cols, right_on=cols)

    return sums


def pixelSize(nside):
    """
    Method to retuen the pixel size

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


class Plot_nsn_vs:
    def __init__(self, data, norm_factor, nside=64):
        """
        class to plot ns vs z or season or ...

        Parameters
        ----------
        data : pandas df
            Data to plot.
        norm_factor : float
            Normalization factor.
        nside: int, optional
            nside healpix parameter. The default is 64.

        Returns
        -------
        None.

        """

        self.data = data
        self.norm_factor = norm_factor
        self.nside = nside

    def plot_versus(self, data, xvar='season', xleg='season',
                    yvar='nsn', yleg='$N_{SN}$', fig=None, ax=None,
                    figTitle='', label=None, xlim=[1, 10],
                    ls='solid', cumul=False, color='k', marker='o'):

        if ax is None:
            fig, ax = plt.subplots(figsize=(14, 9))

        fig.suptitle(figTitle)

        data = data.sort_values(by=[xvar])
        datab = data[yvar]
        if cumul:
            datab = np.cumsum(datab)
        ax.plot(data[xvar], datab, label=label,
                linestyle=ls, marker=marker, color=color, mfc='None', lw=3)
        ax.grid()
        if xlim is not None:
            ax.set_xlim(xlim)

    def plot_nsn_mollview(self, what='season', dbName=''):
        """
        Method to plot the number of SN in Mollweid view

        Parameters
        ----------
        what : TYPE, optional
            DESCRIPTION. The default is 'season'.
        dbName: str, optional
          dbName to display. The default is ''

        Returns
        -------
        None.

        """

        years = self.data[what].unique()

        saveName = '{}_moll'.format(dbName)
        self.Mollview_sum(self.data, addleg='{}'.format(
            dbName), saveName=saveName)

        """
        for year in years:
            idx = self.data[what] == year
            sel = self.data[idx]

            saveName = '{}_moll_{}'.format(dbName, year)
            self.Mollview_sum(sel, addleg='{} \n {} {}'.format(dbName, what, int(year)),
                              saveName=saveName)
        """
        # plt.show()

    def Mollview_sum(self, data, var='nsn',
                     legvar='N$_{SN}$', addleg='', saveName=''):
        """
        Method to plot a Mollweid view for the sum of a variable
        Parameters
        ----------
        data : pandas df
            Data to plot.
        var : str, optional
            Variable to display. The default is 'nsn'.
        legvar : str, optional
            plot legend. The default is 'N$_{SN}$'.
        addleg : str, optional
            Additionnal info for legend. The default is ''.
        saveName : str, optional
            name for the jpeg file. The default is ''.

        Returns
        -------
        None.

        """

        sums = data.groupby(['healpixID']).size().to_frame('nsn').reset_index()
        sums['nsn'] /= self.norm_factor
        print(sums)

        xmin = xmax = np.min(sums[var])
        xmin = 0.1
        xmax = xmax = np.max(sums[var])
        plotMollview(sums, var, legvar, addleg, np.sum,
                     xmin=xmin, xmax=xmax,
                     nside=self.nside, saveName=saveName)


def plot_effi(ax, selb, xvar, yvar, yvar_cut,
              smoothIt, marks, listy, timescale, timeslot, norm_factor):
    """
    Function to plot effi vs z

    Parameters
    ----------
    ax : matplotlib axis
        plot axis.
    selb : pandas df
        Data to plot.
    xvar : str
        x-axis var.
    yvar : str
        y-axis var.
    yvar_cut : float
        y var selection cut.
    smoothIt : bool
        To smooth (spline) displayed curves.
    marks : dict
        Markers used for the plot.
    listy : dict
        linestyle used for the plot.
    timescale : str
        time scale to use (season/year).
    timeslot : int
        time slot for display.
    norm_factor : int
        Normalization factor.

    Returns
    -------
    None.

    """

    df = bin_it_effi(selb, xvar=xvar, yvar=yvar, yvar_cut=yvar_cut,
                     bins=np.arange(0.01, 1.12, 0.1))

    print(df)
    # ax.errorbar(df['z'], df['sigma_mu'], yerr=df['sigma_mu_std'])

    if smoothIt:
        from scipy.interpolate import make_interp_spline
        xnew = np.linspace(
            np.min(df[xvar]), np.max(df[xvar]), 100)
        spl = make_interp_spline(
            df[xvar], df['effi'], k=3)  # type: BSpline
        spl_smooth = spl(xnew)

        ax.plot(xnew, spl_smooth, color='k',
                marker=marks[timeslot], ls=listy[timeslot],
                mfc='None', ms=10, markevery=5,
                label='{} {}'.format(timescale, timeslot))

    else:
        ax.plot(df[xvar], df['effi'], color='k',
                marker=marks[timeslot], ls=listy[timeslot],
                mfc='None', ms=10, markevery=5,
                label='{} {}'.format(timescale, timeslot))
    xmin = 0.2
    xmax = 1.08
    ymin = 0.0
    ymax = None
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    print(df)
    ax.plot([xmin, xmax], [0.95]*2, ls='dashed', color='r')
    ttext = '0.95'
    ax.text(0.3, 0.92, ttext, color='r', fontsize=10)


def plotMollview(data, varName, leg, addleg, op, xmin, xmax,
                 nside=128, saveName=''):
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
    op: operator
      operator to apply to the pixelize data(median, sum, ...)
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
    resleg = op(data[varName])
    if 'nsn' in varName:
        resleg = int(resleg)
    else:
        resleg = np.round(resleg, 2)
    title = '{}: {}'.format(leg, resleg)
    if addleg != '':
        title = '{} - {}'.format(addleg, title)

    hp.mollview(hpxmap, fig=fig, min=xmin, max=xmax, cmap=cmap,
                title=title, nest=True, norm=norm)
    hp.graticule()

    # save plot here
    name = leg.replace(' - ', '_')
    name = name.replace(' ', '_')

    if saveName != '':
        plt.savefig('Plots_pixels/{}.png'.format(saveName))


def plot_nsn_new(data, norm_factor, config, nside,
                 sigma_mu=0.12, timescale='year'):

    # total number of sn per OS/field/timescale

    suma = get_sums_nsn(data, norm_factor, nside, cols=[
        timescale, 'dbName', 'field'])

    # total number of sn per OS/timescale
    sumb = get_sums_nsn(data, norm_factor, nside, cols=[timescale, 'dbName'])

    # total number of sn per OS/timescale
    sumc = get_sums_nsn(data, norm_factor, nside, cols=['dbName'])

    print(sumc)


def get_nsn(data, norm_factor, nside, cols=['year', 'dbName', 'field']):
    """
     Function to get the number of sn and observed area

     Parameters
     ----------
     data : pandas df
         Data to process.
     norm_factor : float
         normalization factor.
     nside : int
         nside healpix parameter.
     cols : list(str), optional
         List of columns to estimate nsn. The default is ['year','dbName','field'].

     Returns
     -------
     pandas df
         result.
     """

    # total number of sn per OS/field/timescale

    sum = get_sums_nsn(data, norm_factor, nside, cols=cols)

    return sum


def plot_versus(data, xvar='season', xleg='season',
                yvar='nsn', yleg='$N_{SN}$', fig=None, ax=None,
                figTitle='', label=None, xlim=[1, 10],
                ls='solid', cumul=False, color='k', marker='o'):

    if ax is None:
        fig, ax = plt.subplots(figsize=(14, 9))

    fig.suptitle(figTitle)

    data = data.sort_values(by=[xvar])
    datab = data[yvar]
    if cumul:
        datab = np.cumsum(datab)
    ax.plot(data[xvar], datab, label=label,
            linestyle=ls, marker=marker, color=color, mfc='None', lw=3)
    ax.grid()
    if xlim is not None:
        ax.set_xlim(xlim)
