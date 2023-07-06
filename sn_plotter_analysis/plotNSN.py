import pandas as pd
import operator
from sn_tools.sn_io import checkDir
from . import plt
import numpy as np
from sn_analysis.sn_calc_plot import select


def plot_NSN(df, xvar='z', xlabel='z', yvar='N', ylabel='NSN',
             bins=np.arange(0.01, 1.1, 0.02),
             norm_factor=1, fig=None, axis=None,
             loopvar='field', dict_sel={}, cumul=False, title='',
             plotDir='', plotName=''):
    """
    Function to plot NSN vs z

    Parameters
    ----------
    df : pandas df
        data to process.
    xvar : str, optional
        x-axis var. The default is 'z'.
    xlabel : str, optional
        x-axis label. The default is 'z'.
    yvar : str, optional
        y-axis var. The default is 'N'.
    ylabel : str, optional
        y-axis label. The default is 'NSN'.
    bins : list(float), optional
        List of bins to draxw a binned-x plot. 
        The default is np.arange(0.01, 1.1, 0.02).
    norm_factor : float, optional
        Normalization factor. The default is 1.
    fig : matplotlib figure, optional
        Figure plot. The default is None.
    axis : matplotlib axis, optional
        axis for the plot. The default is None.
    loopvar : str, optional
        Loop var for the plot. The default is 'field'.
    dict_sel : dict, optional
        Selection dict for the plot. The default is {}.
    cumul : bool, optional
        To plot a cumulative y-var. The default is False.

    Returns
    -------
    None.

    """

    if dict_sel:
        df = select(df, dict_sel['select'])
    fields = df[loopvar].unique()

    if fig is None:
        fig, ax = plt.subplots(figsize=(12, 9))

    fig.subplots_adjust(right=0.8)
    if title != '':
        fig.suptitle(title)

    for field in fields:
        print('field', field)
        idx = df[loopvar] == field
        sel = df[idx]
        if bins is not None:
            sel = bin_it(sel, 'z', bins, norm_factor)
            norm_factor = 1
        yplot = sel[yvar]/norm_factor
        if cumul:
            yplot = np.cumsum(yplot)
        ls = sel['ls'].unique()[0]
        marker = sel['marker'].unique()[0]
        color = sel['color'].unique()[0]
        ax.plot(sel[xvar], yplot, linestyle=ls, marker=marker,
                label=field, mfc='None', color=color)

    ax.grid()
    ax.set_xlabel(xlabel, fontweight='bold')
    ax.set_ylabel(ylabel, fontweight='bold')
    xmin = df[xvar].min()
    xmax = df[xvar].max()
    ax.set_xlim(xmin, xmax)
    # ax.legend(fontsize=12)
    ax.legend(loc='upper center',
              bbox_to_anchor=(1.15, 0.7),
              ncol=1, fontsize=12, frameon=False)
    if plotName != '':
        plt.savefig('{}/{}'.format(plotDir, plotName))


class plotNSN:
    def __init__(self, listDDF, dd, selconfig, selconfig_ref,
                 plotDir='Plots_20230607', fDir=''):
        """
        Class to plot the number of SN

        Parameters
        ----------
        listDDF : str
            List of DDFs to consider.
        dd : pandas df
            list of db to plot (plus plot params).
        selconfig : str
            selection type.
        selconfig_ref : str
            ref for selection.
        plotDir : str, optional
            Plot output dir. The default is 'Plots_20230607'.
        fDir : str, optional
            File output dir. The default is ''.

        Returns
        -------
        None.

        """

        sn_field = pd.read_hdf('{}/sn_field.hdf5'.format(fDir))
        sn_field_season = pd.read_hdf('{}/sn_field_season.hdf5'.format(fDir))

        idx = sn_field_season['field'].isin(listDDF.split(','))

        sn_field_season = sn_field_season[idx]
        sn_tot = sn_field.groupby(['dbName', 'selconfig'])[
            'NSN'].sum().reset_index()
        sn_tot_season = sn_field_season.groupby(['dbName', 'selconfig', 'season'])[
            'NSN'].sum().reset_index()

        self.sn_tot = sn_tot.merge(dd, left_on=['dbName'], right_on=['dbName'])

        self.sn_tot_season = sn_tot_season.merge(
            dd, left_on=['dbName'], right_on=['dbName'])

        self.dict_sel = {}
        self.dict_sel['select'] = [('selconfig', operator.eq, selconfig)]
        self.plotDir = plotDir
        self.tit = '{} \n {}'.format(listDDF, selconfig)
        self.llddf = '_'.join(listDDF.split(','))
        self.selconfig = selconfig
        self.selconfig_ref = selconfig_ref

        checkDir(self.plotDir)

    def plot_NSN_season(self):
        """
        Method to plot NSN vs season

        Returns
        -------
        None.

        """

        plotName = 'nsn_season_{}_{}.png'.format(self.llddf, self.selconfig)
        plot_NSN(self.sn_tot_season,
                 xvar='season', xlabel='season',
                 yvar='NSN', ylabel='$N_{SN}$',
                 bins=None, norm_factor=1, loopvar='dbName',
                 dict_sel=self.dict_sel,
                 title=self.tit, plotDir=self.plotDir, plotName=plotName)

    def plot_NSN_season_cumul(self):
        """
        Method to plot cumulative of Nseasons

        Returns
        -------
        None.

        """
        plotName = 'nsn_sum_season_{}_{}.png'.format(
            self.llddf, self.selconfig)
        plot_NSN(self.sn_tot_season,
                 xvar='season', xlabel='season',
                 yvar='NSN', ylabel='$\Sigma N_{SN}$',
                 bins=None, norm_factor=1,
                 loopvar='dbName', dict_sel=self.dict_sel, cumul=True,
                 title=self.tit, plotDir=self.plotDir, plotName=plotName)

    def plot_observing_efficiency(self):
        """
        Method to plot observing efficiency

        Returns
        -------
        None.

        """

        # plot observing efficiency
        idxa = self.sn_tot_season['selconfig'] == self.selconfig_ref
        sela = self.sn_tot_season[idxa]
        idxb = self.sn_tot_season['selconfig'] == self.selconfig
        selb = self.sn_tot_season[idxb]

        selc = selb.merge(sela, left_on=['season'],
                          right_on=['season'], suffixes=['', '_nosel'])

        print(selc.columns)

        selc['NSN'] /= selc['NSN_nosel']
        selc = selc.sort_values(by=['season'])
        plotName = 'effi_season_{}_{}.png'.format(self.llddf, self.selconfig)
        plot_NSN(selc,
                 xvar='season', xlabel='season',
                 yvar='NSN', ylabel='Observing Efficiency',
                 bins=None, norm_factor=1, loopvar='dbName',
                 dict_sel=self.dict_sel, title=self.tit, plotDir=self.plotDir,
                 plotName=plotName)


def plot_effi(effival, xvar='z', leg='', fig=None, ax=None):
    """
    Function to plot efficiency vs xvar

    Parameters
    ----------
    effival : array
        data to plot.
    xvar : str, optional
        x-axis var. The default is 'z'.
    leg : str, optional
        leg for the plot. The default is ''.
    fig : matplotlib figure, optional
        figure for the plot. The default is None.
    ax : matplotlib axis, optional
        axis for the plot. The default is None.

    Returns
    -------
    None.

    """

    if fig is None:
        fig, ax = plt.subplots(figsize=(10, 8))

    #effival = effi(resa, resb, xvar=xvar, bins=bins)

    x = effival[xvar]
    y = effival['effi']
    yerr = effival['effi_err']
    ax.errorbar(x, y, yerr=yerr, label=leg, marker='o', color='k')
