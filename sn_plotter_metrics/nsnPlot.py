import numpy as np
import pandas as pd
from sn_tools.sn_io import loopStack
import glob
import healpy as hp
import operator

from . import plt
import matplotlib
import matplotlib.patches as mpatches
import tkinter as tk
from tkinter import font as tkFont
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
matplotlib.use('tkagg')


def plot_DDArea(metricValues, forPlot, sntype='faint'):
    """
    Method to plot the 'useful' area wrt dithering

    Parameters
    ----------
    metricValues: pandas df
     data to process
    forPlot: pandas df
      which OS to plot
    sn_type: str
      type of SN to consider (default: faint)

    """

    data = pd.DataFrame(np.copy(metricValues))

    summary_area = data.groupby(['cadence', 'fieldname', 'season']).apply(
        lambda x: useful_area(x)).reset_index()

    # summary_area = summary_area.groupby(['cadence', 'fieldname']).apply(lambda x : pd.DataFrame({'frac_area':[x['frac_area'].median()]})).reset_index()

    summary_area['frac_area'] = summary_area.groupby(
        ['cadence', 'fieldname']).frac_area.transform('median')

    summary_area = summary_area.sort_values(by=['fieldname'])
    print(summary_area)

    fig, ax = plt.subplots()

    varx = 'fieldname'
    vary = 'frac_area'
    for group in np.unique(forPlot['group']):
        idx = forPlot['group'] == group
        sel = forPlot[idx]
        # print(group, sel['dbName'])
        marker = sel['marker'].unique()[0]
        color = sel['color'].unique()[0]

        print('ici', sel['dbName'].str.strip(), summary_area['cadence'])
        selcad = summary_area[summary_area['cadence'].str.strip().isin(
            sel['dbName'].str.strip())]

        selcad = selcad.sort_values(by=['fieldname'], ascending=True)

        # plot
        ax.plot(selcad[varx], selcad[vary], color=color,
                marker=marker, label=group)

    ax.grid()
    ax.set_ylabel('frac_area')
    ax.legend()


def useful_area(grp):
    """
    Method to estimate the 'useful' area wrt dithering

    Parameters
    ----------
    grp: pandas df
     data to process

    Returns
    -------
    panda df with the 'useful' fraction

    """
    npixels = len(np.unique(grp['healpixID']))
    idx = grp['zlim_faint'] > 0.
    sel = grp[idx]

    npixels_useful = len(np.unique(sel['healpixID']))

    return pd.DataFrame({'frac_area': [npixels_useful/npixels]})


def plot_DDSummary(metricValues, forPlot, sntype='faint', fieldNames=['COSMOS'], nside=128):
    """
    Plot to display NSN results for DD fields

    Parameters
    ----------------
    metricValues: numpy array
     array of data to display:
      healpixID: healpixID of the pixel SN
      season: season to plot
      pixRA: RA of the pixel SN
      pixDec: Dec of the pixel SN
      zlim_faint: redshift corresponding to faintest SN (x1=-2.0,color=0.2)
      zlim_medium: redshift corresponding to medium SN (x1=0.0,color=0.0)
      nsn_med_zfaint: number of medium SN with z<zlim_faint
      nsn_med_zmedium: number of medium SN with z<zlim_medium
      nsn_zfaint: number of SN with z<zlim_faint
      nsn_zmedium: number of medium SN with z<zlim_medium
      fieldname: name of the field
      fieldnum: number of the field
      cadence: cadence name
      nside: nside value for healpix tessallation
      pixArea: pixel area
    forPlot: numpy array
     array with a set of plotting information (marker, color) for each cadence:
     dbName: cadence name
     newName: new cadence name
     group: name of the group the cadence is bzelonging to
     Namepl: name of the cadence for plotting
     color: marker color
     marker: marker type
    sntype: str,opt
      type of the supernova (faint or medium) (default: faint) for display
    fieldNames: list(str),opt
      list of fields to consider (default: ['COSMOS'])
    nside: int,opt
      healpix nside value (default: 128)

    Returns
    -----------
    Plot (NSN, zlim)


    """

    # print('selec', sel[['zlim_faint', 'nsn_med_faint']])
    # estimate some stats to display

    data = pd.DataFrame(np.copy(metricValues))
    idx = data['fieldname'].isin(fieldNames)
    idx &= data['zlim_faint'] > 0.
    data = data[idx]

    """
    for season in data['season'].unique():
        idx = data['season'] == season
        sel = data[idx]
        #print(season, sel[['healpixID', 'fieldname', 'zlim_faint']])
        print(season, np.median(sel['zlim_faint']))
        for fieldname in sel['fieldname'].unique():
            idx = sel['fieldname'] == fieldname
            selb = sel[idx]
            print(fieldname, np.median(selb['zlim_faint']))
            #plt.hist(selb['zlim_faint'], histtype='step')
            # plt.show()
    # print(test)
    """
    # plotArea(data, nside)

    """
    summary = data.groupby(['cadence']).agg({'nsn_med_faint': 'sum',
                                             'nsn_med_medium': 'sum',
                                             'zlim_faint': 'median',
                                             'zlim_medium': 'median', }).reset_index()
    """
    summary_fields = data.groupby(['cadence', 'fieldname']).agg({'nsn_med_faint': 'sum',
                                                                 'nsn_med_medium': 'sum',
                                                                 'zlim_faint': 'median',
                                                                 'zlim_medium': 'median', }).reset_index()

    """
    summary_fields_seasons = data.groupby(['cadence', 'fieldname', 'season']).agg({'nsn_med_faint': 'sum',
                                                                                   'nsn_med_medium': 'sum',
                                                                                   'zlim_faint': 'median',
                                                                                   'zlim_medium': 'median', }).reset_index()
    """
    summary_fields_seasons = data.groupby(['cadence', 'fieldname', 'season']).apply(
        lambda x: stat_season(x)).reset_index()

    # print(summary_fields_seasons)
    # print(summary_fields_seasons['zlim_faint_med'].median())

    """
    corresp = dict(zip(['zlim_faint_med', 'zlim_medium_med', 'zlim_faint_weighted', 'zlim_medium_weighted'],
                       ['zlim_medium_med', 'zlim_faint_med', 'zlim_medium_weighted', 'zlim_faint_weighted', 'zlim_medium_weighted']))

    summary = summary_fields_seasons.groupby(['cadence']).apply(
        lambda x: stat_season(x, corresp)).reset_index()
    """

    summary = summary_fields_seasons.groupby(['cadence']).agg({'nsn_med_faint': 'sum',
                                                               'nsn_med_medium': 'sum',
                                                               'zlim_faint_med': 'median',
                                                               'zlim_medium_med': 'median',
                                                               'zlim_faint_weighted': 'median',
                                                               'zlim_medium_weighted': 'median',
                                                               'rms_zlim_faint': 'median',
                                                               'rms_zlim_medium': 'median',
                                                               'rms_zlim_faint_rel': 'median',
                                                               'rms_zlim_medium_rel': 'median'}).reset_index()

    # print(summary)
    """
    summary = data.groupby(['cadence']).agg({'nsn_med_faint': 'sum',
                                             'nsn_med_medium': 'sum',
                                             'zlim_faint': 'median',
                                             'zlim_medium': 'median', }).reset_index()
    """
    # print(summary_fields_seasons)
    # print('aiaiai',summary.columns)
    # change some of the type for printing
    summary.round({'zlim_faint_med': 2, 'zlim_medium_med': 2})
    summary['nsn_med_faint'] = summary['nsn_med_faint'].astype(int)
    summary['nsn_med_medium'] = summary['nsn_med_medium'].astype(int)

    # plot the results

    # per field and per season
    # plotNSN(summary_fields_seasons, forPlot, sntype='faint', ztype='med')
    # plotNSN(summary_fields_seasons, forPlot, sntype='faint', ztype='weighted')

    # plotNSN(data, forPlot, sntype=sntype, varx='zlim_faint')

    # per field, for all seasons
    # plotNSN(summary_fields, forPlot, sntype=sntype, varx='zlim_faint')
    # Summary plot: one (NSN,zlim) per cadence (sum for NSN, median zlim over the fields/seasons)

    #plotNSN(summary, forPlot, sntype=sntype, varx='zlim_faint_med')
    """
    for fieldname in summary_fields['fieldname'].unique():
        idx = summary_fields['fieldname'] == fieldname
        sel = summary_fields[idx]
        plotNSN(sel, forPlot, varx='zlim_faint')
    """
    plotNSN(summary_fields, forPlot, varx='zlim_faint')

    # plotNSN(summary, forPlot, sntype=sntype, varx='zlim_faint_weighted')

    # this is for dithering plots wrt wRMS
    """
    idx = summary['cadence'] == 'ddf_dither0.00_v1.7_10yrs'
    norm = summary[idx]
    plotNSN(summary, forPlot,
            varx='rms_zlim_{}'.format(sntype),
            vary='nsn_med_{}'.format(sntype),
            legx='weighted RMS($z_{lim}$)',
            legy='$N_{SN}/N_{SN}^{no dither} (z<)$',
            norm=norm['nsn_med_{}'.format(sntype)].item())

    plotNSN(summary, forPlot,
            varx='rms_zlim_{}'.format(sntype),
            vary='zlim_{}_med'.format(sntype),
            legx='weighted RMS($z_{lim}$)',
            legy='$\Delta z_{lim}=z_{lim}^{no dither}-z_{lim}$',
            norm=norm['zlim_{}_med'.format(sntype)].item(), op=operator.sub)
    """


def stat_season(grp,
                corresp=dict(zip(['zlim_faint_med', 'zlim_medium_med', 'zlim_faint_weighted', 'zlim_medium_weighted'],
                                 ['zlim_faint', 'zlim_medium', 'zlim_faint', 'zlim_medium']))):
    """
    Method to estimate weighted mean and rms of zlim

    Parameters
    ----------------
    grp: pandas df
      data to process
    corresp: dict
      matching dict for zlim vars

    """
    dictres = {}

    for vv in ['faint', 'medium']:

        nsn = 'nsn_med_{}'.format(vv)
        zlim_med = 'zlim_{}_med'.format(vv)
        zlim_weighted = 'zlim_{}_weighted'.format(vv)
        rms_zlim = 'rms_zlim_{}'.format(vv)
        rms_zlim_rel = 'rms_zlim_{}_rel'.format(vv)
        weights = grp[nsn]
        dictres[nsn] = [grp[nsn].sum()]
        dictres[zlim_med] = [grp[corresp[zlim_med]].median()]
        zlim_weight = np.sum(weights*grp[corresp[zlim_weighted]])/weights.sum()
        dictres[zlim_weighted] = [zlim_weight]
        std = np.sum(weights*(grp[corresp[zlim_weighted]] -
                              zlim_weight)**2)/np.sum(weights)
        dictres[rms_zlim] = [np.sqrt(std)]
        dictres[rms_zlim_rel] = [np.sqrt(std)/zlim_weight]

    return pd.DataFrame(dictres)


def plotNSN(summary, forPlot,
            varx='zlim_faint_med',
            vary='nsn_med_faint',
            legx='$z_{faint}$',
            legy='$N_{SN} (z<)$',
            norm=1,
            op=operator.truediv):
    """
    Plot NSN vs redshift limit

    Parameters
    ----------------
    summary: pandas Dataframe
     data to display:
      cadence: name of the cadence
      zlim_faint: redshift corresponding to faintest SN (x1=-2.0,color=0.2)
      zlim_medium: redshift corresponding to medium SN (x1=0.0,color=0.0)
      nsn_zfaint: number of SN with z<zlim_faint
      nsn_zmedium: number of medium SN with z<zlim_medium
    forPlot: numpy array
      array with a set of plotting information (marker, color) for each cadence:
      dbName: cadence name
      newName: new cadence name
      group: name of the group the cadence is bzelonging to
      Namepl: name of the cadence for plotting
      color: marker color
      marker: marker type
    sntype: str,opt
      type of the supernova (faint or medium) for the display (default: faint)

    Returns
    -----------
    Plot (NSN, zlim)


    """

    fontsize = 15
    fig, ax = plt.subplots()
    # varx = 'zlim_{}_{}'.format(sntype, ztype)
    #vary = 'nsn_med_{}'.format(sntype)
    xshift = 1.0
    yshift = 1.01

    for group in np.unique(forPlot['group']):
        idx = forPlot['group'] == group
        sel = forPlot[idx]
        # print(group, sel['dbName'])
        marker = sel['marker'].unique()[0]
        color = sel['color'].unique()[0]

        # print('ici', sel['dbName'].str.strip(), summary['cadence'])
        selcad = summary[summary['cadence'].str.strip().isin(
            sel['dbName'].str.strip())]

        # plot
        ax.plot(selcad[varx], op(selcad[vary], norm), color=color,
                marker=marker, lineStyle='None')

        # get the centroid of the data and write it
        centroid_x = selcad[varx].mean()
        centroid_y = op(selcad[vary], norm).mean()
        ax.text(xshift*0.99*centroid_x, yshift *
                1.01*centroid_y, group, color=color, fontsize=fontsize-3)

    ax.grid()
    ax.set_xlabel(legx)
    ax.set_ylabel(legy)

    """
    fig.text(0.8, 0.8, 'Preliminary',
             fontsize=25, color='blue',
             ha='right', va='bottom', alpha=0.5)
    """


def mscatter(x, y, ax=None, m=None, **kw):
    import matplotlib.markers as mmarkers
    ax = ax or plt.gca()
    sc = ax.scatter(x, y, **kw)
    if (m is not None) and (len(m) == len(x)):
        paths = []
        for marker in m:
            if isinstance(marker, mmarkers.MarkerStyle):
                marker_obj = marker
            else:
                marker_obj = mmarkers.MarkerStyle(marker)
            path = marker_obj.get_path().transformed(
                marker_obj.get_transform())
            paths.append(path)
        sc.set_paths(paths)
    return sc


class NSNAnalysis:
    def __init__(self, dbInfo,
                 metricName='NSN', fieldType='WFD',
                 nside=64,
                 x1=-2.0, color=0.2, npixels=-1):
        """
        class to analyze results from NSN metric

        Parameters
        ---------------
        dbInfo: pandas df
          info from observing strategy (simutype, simuname,dbDir,dbName, , color, marker)
        metricName: str
          name of the metric used to generate the files
        fieldType: str, opt
          field type (DD, WFD) (default: WFD)
        nside: int,opt
          healpix nside parameter (default: 64)
        x1: float, opt
          x1 SN value (default: -2.0)
        color: float, opt
          color SN value (default: 0.2)
        npixels: int, opt
          total number of pixels for this strategy

        """

        self.nside = nside
        self.npixels = npixels
        self.pixel_area = hp.nside2pixarea(nside, degrees=True)

        self.sntype = 'faint'
        if x1 == 0. and color == 0.:
            self.sntype = 'medium'

        self.dbInfo = dbInfo

        # loading data (metric values)
        search_path = '{}/{}/{}/*NSNMetric_{}*_nside_{}_*.hdf5'.format(
            dbInfo['dirFile'], dbInfo['dbName'], metricName, fieldType, nside)
        print('looking for', search_path)
        fileNames = glob.glob(search_path)
        # fileName='{}/{}_CadenceMetric_{}.npy'.format(dirFile,dbName,band)
        print(fileNames)
        if len(fileNames) > 0:
            self.data_summary = self.process(fileNames)
        else:
            print('Missing files for', dbInfo['dbName'])
            self.data_summary = None

    def process(self, fileNames):
        """
        Method to process metric values from files

        Parameters
        ---------------
        fileNames: list(str)
          list of files to process

        Returns
        ----------
        resdf: pandas df with a summary of metric infos

        """
        metricValues = np.array(loopStack(fileNames, 'astropyTable'))

        """
        df = pd.DataFrame(np.copy(metricValues))
        df = df.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)
        # self.data = df.loc[:,~df.columns.str.contains('mask', case=False)]
       """

        idx = metricValues['status_{}'.format(self.sntype)] == 1
        # idx &= metricValues['healpixID'] >= 48000
        # idx &= metricValues['healpixID'] <= 49000

        idx &= metricValues['zlim_{}'.format(self.sntype)] > 0.
        idx &= metricValues['nsn_med_{}'.format(self.sntype)] > 0.

        # self.plot_season(metricValues[idx], varName='nsn_med')

        self.data = pd.DataFrame(metricValues[idx])
        self.data = self.data.applymap(
            lambda x: x.decode() if isinstance(x, bytes) else x)
        print('data', self.data[['healpixID', 'pixRA', 'pixDec', 'zlim_{}'.format(self.sntype),
                                 'nsn_med_{}'.format(self.sntype), 'nsn', 'season']], self.data.columns)
        print(len(np.unique(self.data[['healpixID', 'season']])))
        self.ratiopixels = 1
        self.npixels_eff = len(self.data['healpixID'].unique())
        if self.npixels > 0:
            self.ratiopixels = float(
                npixels)/float(self.npixels_eff)

        # zlim = self.zlim_med()
        nsn, sig_nsn = self.nSN_tot()
        nsn_extrapol = int(np.round(nsn*self.ratiopixels))

        """
        test = self.data.apply(lambda x: self.ana_filters(
            x),axis=1,result_type='expand')
        for b in 'grizy':
            tt=[]
            for vv in self.data.columns:
                if 'N_' in vv and b in vv:
                    tt.append(vv)
            self.data['N_{}_tot'.format(b)] = self.data[tt].sum(axis=1)
            print(self.data[tt])
            print('sum',self.data['N_{}_tot'.format(b)])

        """
        bandstat = ['u', 'g', 'r', 'i', 'z', 'y', 'gr', 'gi',
                    'gz', 'iz', 'uu', 'gg', 'rr', 'ii', 'zz', 'yy']

        for b in bandstat:
            self.data['cadence_{}'.format(
                b)] = self.data['season_length']/self.data['N_{}'.format(b)]

        meds = self.data.groupby(['healpixID']).median().reset_index()
        meds = meds.round({'zlim_{}'.format(self.sntype): 5})
        med_meds = meds.median()

        resdf = pd.DataFrame(
            [med_meds['zlim_{}'.format(self.sntype)]], columns=['zlim'])
        resdf['nsn'] = [nsn]
        resdf['sig_nsn'] = [sig_nsn]
        resdf['nsn_extra'] = [nsn_extrapol]
        resdf['dbName'] = self.dbInfo['dbName']
        resdf['simuType'] = self.dbInfo['simuType']
        resdf['simuNum'] = self.dbInfo['simuNum']
        resdf['family'] = self.dbInfo['family']
        resdf['color'] = self.dbInfo['color']
        resdf['marker'] = self.dbInfo['marker']
        resdf['cadence'] = [med_meds['cadence']]
        resdf['season_length'] = [med_meds['season_length']]
        resdf['N_total'] = [med_meds['N_total']]
        resdf['survey_area'] = self.npixels_eff*self.pixel_area
        resdf['nsn_per_sqdeg'] = resdf['nsn']/resdf['survey_area']

        for b in bandstat:
            resdf['N_{}'.format(b)] = [med_meds['N_{}'.format(b)]]
            resdf['cadence_{}'.format(b)] = [med_meds['cadence_{}'.format(b)]]
            """
            resdf['N_{}_tot'.format(ba)] = [med_meds['N_{}_tot'.format(ba)]]
            resdf['N_{}'.format(ba)] = [
                                med_meds['N_{}'.format(ba)]/resdf['N_total']]
            for bb in 'grizy':
                combi=''.join(sorted('{}{}'.format(ba,bb)))
                # resdf['cadence_{}{}'.format(ba,bb)] = [med_meds['cadence_{}{}'.format(ba,bb)]]
                resdf['N_{}'.format(combi)] = [
                                    med_meds['N_{}'.format(combi)]/resdf['N_total']]
                for bc in 'grizy':
                    combi=''.join(sorted('{}{}{}'.format(ba,bb,bc)))
                    # resdf['cadence_{}{}{}'.format(ba,bb,bc)] = [med_meds['cadence_{}{}{}'.format(ba,bb,bc)]]
                    resdf['N_{}'.format(combi)] = [
                                        med_meds['N_{}'.format(combi)]/resdf['N_total']]
            """

        return resdf

    def ana_filters(self, grp):

        print('io', grp['N_filters_night'])
        filter_night = grp['N_filters_night'].split('/')
        print(filter_night)
        r = []
        for vv in filter_night:
            if vv != '':
                vs = vv.split('*')
                print('what', vs)
                r.append((int(vs[0]), vs[1]))

        print(r)
        res = pd.DataFrame(r, columns=['Nvisits', 'filter'])

        print(res['Nvisits'].sum(), grp['N_total'])

        return pd.DataFrame({'test': grp})

    def plot_season(self, metricValues, varName='status', op=np.sum):
        """
        Method to plot seasonal results

        Parameters
        ---------------
        metricValues: pandas df
           data to plot
        varName: str, opt
            variable to plot(default: status)

        """

        # loop on seasons and plot (Mollview) the variable mean per pixel
        for season in np.unique(metricValues['season']):
            idx = metricValues['season'] == season
            sel = metricValues[idx]
            leg = '{} - season {}'.format(varName, int(season))
            self.plotMollview(sel, varName, leg, op,
                              np.min(sel[varName]), np.max(sel[varName]))

        # now for all seasons: plot the max of the variable per pixel
        leg = '{} - all seasons'.format(varName)
        """
        go = pd.DataFrame(metricValues).groupby(
            ['healpixID']).max().reset_index()
        for io, row in go.iterrows():
            print(row[['healpixID', varName]].values)
        """
        self.plotMollview(metricValues, varName, leg, np.median,
                          np.median(metricValues[varName]), np.median(metricValues[varName]))

        plt.show()

    def zlim_med(self):
        """
        Method to estimate the median redshift limit over the pixels

        Returns
        ----------
        median zlim(float)
        """
        meds = self.data.groupby(['healpixID']).median().reset_index()
        meds = meds.round({'zlim': 2})

        return meds['zlim'].median()

    def nSN_tot(self):
        """
        Method to estimate the total number of supernovae(and error)

        Returns
        -----------
        nsn, sig_nsn: int, int
          number of sn and sigma
        """
        sums = self.data.groupby(['healpixID']).sum().reset_index()

        """
        for ii, row in sums.iterrows():
            print(row[['healpixID', 'nsn_med']].values)

        sums['nsn_med'] = sums['nsn_med'].astype(int)
        """
        return sums['nsn_med_{}'.format(self.sntype)].sum(), int(np.sqrt(sums['err_nsn_med_{}'.format(self.sntype)].sum()))

    def Mollview_median(self, var='zlim', legvar='zlimit'):
        """
        Method to plot a Mollweid view for the median of a variable

        Parameters
        --------------
        var: str,opt
          variable to show (default: nsn_med)
        legvar: str, opt
           name for title of the plot (default: NSN)

        """

        # this is to estimate the median zlim over the sky
        meds = self.data.groupby(['healpixID']).median().reset_index()
        meds = meds.round({var: 2})
        self.plotMollview(meds, var, legvar, np.median,
                          xmin=0.000001, xmax=np.max(meds[var]))

    def Mollview_sum(self, var='nsn_med', legvar='NSN'):
        """
        Method to plot a Mollweid view for the sum of a variable

        Parameters
        --------------
        var: str,opt
          variable to show (default: nsn_med)
        legvar: str, opt
           name for title of the plot (default: NSN)

        """

        sums = self.data.groupby(['healpixID']).sum().reset_index()

        self.plotMollview(sums, var, legvar, np.sum,
                          xmin=np.min(sums[var]), xmax=np.max(sums[var]))

    def plot(self):
        """
        Method to plot two Mollview of the metric results:
        - redshift limit
        - number of well-sampled supernovae

        """

        # this is to estimate the median zlim over the sky
        meds = self.data.groupby(['healpixID']).median().reset_index()
        meds = meds.round({'zlim': 2})
        self.plotMollview(meds, 'zlim', 'zlimit', np.median,
                          xmin=0.01, xmax=np.max(meds['zlim'])+0.1)

        # this is to plot the total number of SN (per pixels) over the sky
        sums = self.data.groupby(['healpixID']).sum().reset_index()

        self.plotMollview(sums, 'nsn_med', 'NSN', np.sum,
                          xmin=np.min(sums['nsn_med']), xmax=np.max(sums['nsn_med']))

        """
        self.plotMollview(self.data, 'healpixID', 'healpixID', np.mean,
                          xmin=0.0001, xmax=np.max(self.data['healpixID'])+1)
        """

    def plotMollview(self, data, varName, leg, op, xmin, xmax):
        """
        Method to display results as a Mollweid map

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

        """
        npix = hp.nside2npix(self.nside)

        hpxmap = np.zeros(npix, dtype=np.float)
        hpxmap = np.full(hpxmap.shape, 0.)
        hpxmap[data['healpixID'].astype(
            int)] += data[varName]

        norm = plt.cm.colors.Normalize(xmin, xmax)
        cmap = plt.cm.jet
        cmap.set_under('w')
        resleg = op(data[varName])
        if 'nsn' in varName:
            resleg = int(resleg)
        else:
            resleg = np.round(resleg, 2)
        title = '{}: {}'.format(leg, resleg)

        hp.mollview(hpxmap, min=xmin, max=xmax, cmap=cmap,
                    title=title, nest=True, norm=norm)
        hp.graticule()

        # save plot here
        name = leg.replace(' - ', '_')
        name = name.replace(' ', '_')

        plt.savefig('Plots_pixels/Moll_{}.png'.format(name))

    def plotCorrel(self, datax, datay, x=('', ''), y=('', '')):
        """
        Method for 2D plots

        Parameters
        ---------------
        datax: pandas df
          data used to display - x-axis
        datay: pandas df
          data used to display - y-axis
        x : tuple(str)
          tuple for the x-axis variable
          the first value is the colname in self.data
          the second value is the x-axis corresponding legend
         y : tuple(str)
          tuple for the y-axis variable
          the first value is the colname in self.data
          the second value is the y-axis corresponding legend

        """
        from scipy.stats import pearsonr
        fig, ax = plt.subplots()
        fig.suptitle(self.dbInfo['dbName'])

        varx = datax[x[0]]
        vary = datay[y[0]]
        ax.plot(varx, vary, 'k.')
        corr, _ = pearsonr(varx, vary)
        print('Pearsons correlation: %.3f' % corr)
        cov = np.cov(varx, vary)
        print('Covariance matrix', cov)

        ax.set_xlabel(x[1])
        ax.set_ylabel(y[1])


class PlotSummary_Annot:
    """
    class to plot the metric (nSN,zlim)
    this is an interactive plot (ie infos on scattered points using the mouse pointer)

    Parameters
    ---------------
    resdf: pandas df
      data to plot
    fig: matplotlib figure
     figure where to plot
    ax: matplotlib axes
     axes where to plot
    hlist: list
      list of OS to highlight (will appear in gold)

    """

    def __init__(self, resdf, hlist=[]):

        self.fig, self.ax = plt.subplots(figsize=(12, 8))
        # self.ax = ax
        self.datadf = resdf
        x = self.datadf['zlim'].to_list()
        y = self.datadf['nsn'].to_list()
        c = self.datadf['color'].to_list()
        m = self.datadf['marker'].to_list()
        # self.fig, self.ax = plt.subplots(figsize=(14, 8))

        self.sc = mscatter(x, y, ax=self.ax, c=c, m=m, s=100)

        # this is for the interactive part
        self.annot = self.ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                                      bbox=dict(boxstyle="round", fc="w"),
                                      arrowprops=dict(arrowstyle="->"))
        self.annot.set_visible(False)
        self.ax.grid()
        self.ax.set_xlabel('$z_{faint}$', fontsize=20)
        self.ax.set_ylabel('$N_{SN}(z\leq z_{faint})$', fontsize=20)
        self.ax.tick_params(labelsize=15)
        patches = []
        for col in np.unique(self.datadf['color']):
            idx = self.datadf['color'] == col
            tab = self.datadf[idx]
            lab = '{}_{}'.format(
                np.unique(tab['simuType']).item(), np.unique(tab['simuNum']).item())
            patches.append(mpatches.Patch(color=col, label=lab))
        self.ax.legend(handles=patches, fontsize=15, loc='upper left')

        # highlight here
        if hlist:
            resh = self.select(self.datadf, hlist)
            resh['color'] = 'gold'
            x = resh['zlim'].to_list()
            y = resh['nsn'].to_list()
            c = resh['color'].to_list()
            m = resh['marker'].to_list()
            mscatter(x, y, ax=self.ax, c=c, m=m, s=90)

        self.fig.canvas.mpl_connect("motion_notify_event", self.hover)

        plt.show()

    def update_annot(self, ind):
        """
        Method to update annotation on the plot

        Parameters
        ---------------
        ind: dict
           point infos

        """

        pos = self.sc.get_offsets()[ind["ind"][0]]
        self.annot.xy = pos
        idx = np.abs(self.datadf['zlim']-pos[0]) < 1.e-8
        idx &= np.abs(self.datadf['nsn']-pos[1]) < 1.e-8
        sel = self.datadf[idx]
        color = sel['color'].values.item()

        text = "{} \n family: {}".format(
            sel['dbName'].values.item(), sel['family'].values.item())

        self.annot.set_text(text)
        # self.annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
        self.annot.get_bbox_patch().set_facecolor(color)
        self.annot.get_bbox_patch().set_alpha(0.4)

    def hover(self, event):
        """
        Method to plot annotation when pointer is on a point

        Parameters
        ---------------
        event: matplotlib.backend_bases.MouseEvent
          mouse pointer?
        """
        vis = self.annot.get_visible()
        if event.inaxes == self.ax:
            cont, ind = self.sc.contains(event)
            if cont:
                self.update_annot(ind)
                self.annot.set_visible(True)
                self.fig.canvas.draw_idle()
            else:
                if vis:
                    self.annot.set_visible(False)
                    self.fig.canvas.draw_idle()

    def select(self, resdf, strfilt=['noddf']):
        """
        Function to remove OS according to their names

        Parameters
        ---------------
        resdf: pandas df
        data to process
        strfilt: list(str),opt
         list of strings used to select OS (default: ['noddf']

        """

        for vv in strfilt:
            idx = resdf['dbName'].str.contains(vv)
            resdf = pd.DataFrame(resdf[idx])

        return resdf


class NSN_zlim_GUI:
    """
    class to display the (nsn,zlim) metric as an interactive plot in a GUI

    Parameters
    ---------------
    resdf: pandas df
      data to display

    """

    def __init__(self, resdf):

        self.resdf = resdf

        # build the GUI here
        root = tk.Tk()
        # figure where the plots will be drawn
        self.fig = plt.Figure(figsize=(14, 8), dpi=100)
        self.ax = self.fig.add_subplot(111)
        leg = 'days$^{-1}$'
        # self.fig.suptitle('(nSN,zlim) supernovae metric', fontsize=15)
        self.fig.subplots_adjust(right=0.75)
        # self.ax.set_xlim(self.zmin, self.zmax)
        # define the figure canvas here
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH)
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=False)

        self.toolbar = NavigationToolbar2Tk(self.canvas, root)
        self.toolbar.update()
        # self.ax.cla()

        # plot the metric here
        # PlotSummary_Annot(self.resdf, fig=self.fig, ax=self.ax)
        PlotSummary_Annot(self.fig, self.ax, self.resdf)
        # common font
        helv36 = tkFont.Font(family='Helvetica', size=15, weight='bold')

        # building the GUI
        # frame
        button_frame = tk.Frame(master=root, bg="white")
        button_frame.pack(fill=tk.X, side=tk.BOTTOM, expand=False)
        button_frame.place(relx=.88, rely=.5, anchor="c")

        button_frameb = tk.Frame(master=root, bg="white")
        button_frameb.pack(fill=tk.X, side=tk.BOTTOM, expand=False)
        button_frameb.place(relx=.92, rely=.9, anchor="c")

        # entries
        ents = self.make_entries(button_frame, font=helv36)

        # buttons
        heightb = 3
        widthb = 6

        update_button = tk.Button(
            button_frameb, text="Update", command=(lambda e=ents: self.updateGUI(e)),
            bg='yellow', height=heightb, width=widthb, fg='blue', font=helv36)
        quit_button = tk.Button(button_frameb, text="Quit",
                                command=root.quit, bg='yellow',
                                height=heightb, width=widthb, fg='black', font=helv36)

        button_frame.columnconfigure(0, weight=1)
        button_frame.columnconfigure(1, weight=1)

        update_button.grid(row=1, column=0, sticky=tk.W+tk.E)
        quit_button.grid(row=1, column=1, sticky=tk.W+tk.E)

        root.mainloop()

    def updateGUI(self, entries):
        """
        Method to update the figure according to request made on entries
        The metric (nsn, zlim) will be plotted taking entries value into account.

        Parameters
        ---------------
        entries: dict of tk.Entry
        """

        # reset axes
        self.ax.cla()

        # get the list of OS to drop (from GUI)
        droplist = entries['drop'].get('1.0', tk.END).split('\n')
        droplist = ' '.join(droplist).split()
        resfi = self.filter(self.resdf, droplist)

        # get the list of OS to highligth (from GUI)
        highlightlist = entries['highlight'].get('1.0', tk.END).split('\n')
        highlightlist = ' '.join(highlightlist).split()

        # display here
        # PlotSummary_Annot(resfi, fig=self.fig, ax=self.ax, hlist=highlightlist)
        PlotSummary_Annot(self.fig, self.ax, resfi, hlist)
        # update canvas
        # self.ax.set_xlim(self.zmin, self.zmax)
        self.canvas.draw()

    def make_entries(self, frame, font):
        """
        Method to define entries to the GUI
        Parameters
        ---------------
        frame: tk.Frame
          frame where entries will be located
        font: tk.Font
          font for entries
        Returns
        ----------
        entries: dict
          dict of tk.Entry
        """

        tk.Label(frame, text='drop', bg='white',
                 fg='blue', font=font).grid(row=0)
        tk.Label(frame, text='highlight', bg='white',
                 fg='red', font=font).grid(row=1)

        entries = {}

        # this is to drop OS
        S = tk.Scrollbar(frame)
        entries['drop'] = tk.Text(frame, height=4, width=20, font=font)
        S.grid(row=0, column=1)
        entries['drop'].grid(ipady=80, row=0, column=1)

        S.config(command=entries['drop'].yview)
        entries['drop'].config(yscrollcommand=S.set)
        defa = 'noddf\nfootprint_stuck_rolling\nweather\nwfd_depth_scale'
        entries['drop'].insert(tk.END, defa)

        # this is to highlight OS
        Sh = tk.Scrollbar(frame)
        entries['highlight'] = tk.Text(frame, height=4, width=20, font=font)
        Sh.grid(row=1, column=1)
        entries['highlight'].grid(ipady=20, row=1, column=1)

        Sh.config(command=entries['highlight'].yview)
        entries['highlight'].config(yscrollcommand=S.set)
        entries['highlight'].insert(tk.END, '')

        return entries

    def filter(self, resdf, strfilt=['_noddf']):
        """
        Function to remove OS according to their names

        Parameters
        ---------------
        resdf: pandas df
        data to process
        strfilt: list(str),opt
         list of strings used to remove OS (default: ['_noddf']

        """

        for vv in strfilt:
            idx = resdf['dbName'].str.contains(vv)
            resdf = pd.DataFrame(resdf[~idx])

        return resdf


class NSN_zlim_GUI:
    def __init__(self, resdf):
        self.resdf = resdf

        # build the GUI here
        root = tk.Tk()
        # figure where the plots will be drawn
        self.fig = plt.Figure(figsize=(16, 8), dpi=100)
        self.ax = self.fig.add_subplot(111)
        # leg = 'days$^{-1}$'
        # self.fig.suptitle('(nSN,zlim) supernovae metric', fontsize=15)
        self.fig.subplots_adjust(right=0.75)
        # self.ax.set_xlim(self.zmin, self.zmax)
        # define the figure canvas here
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH)
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=False)

        self.toolbar = NavigationToolbar2Tk(self.canvas, root)
        self.toolbar.update()
        # self.ax.cla()

        # plot the number of visits vs z
        self.plotMetric(self.resdf)

        # common font
        helv36 = tkFont.Font(family='Helvetica', size=15, weight='bold')

        # building the GUI
        # frame
        button_frame = tk.Frame(master=root, bg="white")
        button_frame.pack(fill=tk.X, side=tk.BOTTOM, expand=False)
        button_frame.place(relx=.88, rely=.5, anchor="c")

        button_frameb = tk.Frame(master=root, bg="white")
        button_frameb.pack(fill=tk.X, side=tk.BOTTOM, expand=False)
        button_frameb.place(relx=.92, rely=.9, anchor="c")

        # entries
        ents = self.make_entries(button_frame, font=helv36)

        # buttons
        heightb = 3
        widthb = 6

        update_button = tk.Button(
            button_frameb, text="Update", command=(lambda e=ents: self.updateGUI(e)),
            bg='yellow', height=heightb, width=widthb, fg='blue', font=helv36)
        quit_button = tk.Button(button_frameb, text="Quit",
                                command=root.quit, bg='yellow',
                                height=heightb, width=widthb, fg='black', font=helv36)

        button_frame.columnconfigure(0, weight=1)
        button_frame.columnconfigure(1, weight=1)
        # button_frame.columnconfigure(2, weight=1)

        """
        nvisits_button.grid(row=2, column=0, sticky=tk.W+tk.E)
        z_button.grid(row=2, column=1, sticky=tk.W+tk.E)
        """
        update_button.grid(row=1, column=0, sticky=tk.W+tk.E)
        quit_button.grid(row=1, column=1, sticky=tk.W+tk.E)

        root.mainloop()

    def updateGUI(self, entries):
        """
        Method to update the figure according to request made on entries
        The number of visits will be plotted here
        Parameters
        ---------------
        entries: dict of tk.Entry
        """

        # reset axes
        self.ax.cla()
        # plot Nvisits vs z
        # filter cadences here

        droplist = entries['drop'].get('1.0', tk.END).split('\n')
        droplist = ' '.join(droplist).split()
        resfi = self.filter(self.resdf, droplist)
        highlightlist = entries['highlight'].get('1.0', tk.END).split('\n')
        highlightlist = ' '.join(highlightlist).split()

        norm = entries['norm'].get()
        print('hhh', norm)
        self.plotMetric(resfi, highlightlist, norm)

        # update canvas
        # self.ax.set_xlim(self.zmin, self.zmax)
        self.canvas.draw()

    def make_entries(self, frame, font):
        """
        Method to define entries to the GUI
        Parameters
        ---------------
        frame: tk.Frame
          frame where entries will be located
        font: tk.Font
          font for entries
        Returns
        ----------
        entries: dict
          dict of tk.Entry
        """

        tk.Label(frame, text='drop', bg='white',
                 fg='blue', font=font).grid(row=0)
        tk.Label(frame, text='highlight', bg='white',
                 fg='red', font=font).grid(row=1)
        tk.Label(frame, text='norm', bg='white',
                 fg='black', font=font).grid(row=2)

        entries = {}

        # this is to drop OS
        S = tk.Scrollbar(frame)
        entries['drop'] = tk.Text(frame, height=4, width=20, font=font)
        S.grid(row=0, column=1)
        entries['drop'].grid(ipady=80, row=0, column=1)

        S.config(command=entries['drop'].yview)
        entries['drop'].config(yscrollcommand=S.set)
        defa = 'noddf\nfootprint_stuck_rolling\nweather\nwfd_depth_scale'
        entries['drop'].insert(tk.END, defa)

        # this is to highlight OS
        Sh = tk.Scrollbar(frame)
        entries['highlight'] = tk.Text(frame, height=4, width=20, font=font)
        Sh.grid(row=1, column=1)
        entries['highlight'].grid(ipady=20, row=1, column=1)

        Sh.config(command=entries['highlight'].yview)
        entries['highlight'].config(yscrollcommand=S.set)
        entries['highlight'].insert(tk.END, '')

        entries['norm'] = tk.Entry(frame, width=20, font=font)
        entries['norm'].insert(10, "")
        entries['norm'].grid(row=2, column=1)

        return entries

    def filter(self, resdf, strfilt=['_noddf']):
        """
        Function to remove OS according to their names

        Parameters
        ---------------
        resdf: pandas df
        data to process
        strfilt: list(str),opt
         list of strings used to remove OS (default: ['_noddf']

        """

        for vv in strfilt:
            idx = resdf['dbName'].str.contains(vv)
            resdf = pd.DataFrame(resdf[~idx])

        return resdf

    def select(self, resdf, strfilt=['noddf']):
        """
        Function to remove OS according to their names

        Parameters
        ---------------
        resdf: pandas df
        data to process
        strfilt: list(str),opt
         list of strings used to select OS (default: ['noddf']

        """

        for vv in strfilt:
            idx = resdf['dbName'].str.contains(vv)
            resdf = pd.DataFrame(resdf[idx])

        return resdf

    def plotMetric(self, resdf, hlist=[], norm=''):

        nOS = len(resdf)
        title = '(nSN,zlim) supernovae metric - {} OS'.format(nOS)

        resdf['nsn_norm'] = resdf['nsn']
        # normalize nsn
        ido = resdf['dbName'] == norm
        if len(resdf[ido]) > 0:
            knorm = resdf[ido]['nsn'].values.item()
            resdf['nsn_norm'] = resdf['nsn']/knorm

        self.fig.suptitle(title, fontsize=15)
        self.resdfa = resdf
        x = resdf['zlim'].to_list()
        y = resdf['nsn_norm'].to_list()
        c = resdf['color'].to_list()
        m = resdf['marker'].to_list()
        # self.fig, self.ax = plt.subplots(figsize=(14, 8))

        self.sc = mscatter(x, y, ax=self.ax, c=c, m=m, s=100)
        # self.sc = self.ax.plot(x, y, color=c)
        self.annot = self.ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                                      bbox=dict(boxstyle="round", fc="w"),
                                      arrowprops=dict(arrowstyle="->"))
        self.annot.set_visible(False)

        self.ax.grid()
        self.ax.set_xlabel('$z_{faint}$', fontsize=20)
        self.ax.set_ylabel('$N_{SN}(z\leq z_{faint})$', fontsize=20)
        self.ax.tick_params(labelsize=15)
        patches = []
        for col in np.unique(resdf['color']):
            idx = resdf['color'] == col
            tab = resdf[idx]
            lab = '{}_{}'.format(
                np.unique(tab['simuType']).item(), np.unique(tab['simuNum']).item())
            patches.append(mpatches.Patch(color=col, label=lab))
        self.ax.legend(handles=patches, fontsize=15, loc='upper left')

        if hlist:
            resh = self.select(resdf, hlist)
            resh['color'] = 'gold'
            x = resh['zlim'].to_list()
            y = resh['nsn_norm'].to_list()
            c = resh['color'].to_list()
            m = resh['marker'].to_list()
            mscatter(x, y, ax=self.ax, c=c, m=m, s=90)

        self.fig.canvas.mpl_connect("motion_notify_event", self.hover)

        # plt.show()

    def update_annot(self, ind):
        """
        Method to update annotation on the plot

        Parameters
        ---------------
        ind: dict
           point infos

        """
        pos = self.sc.get_offsets()[ind["ind"][0]]
        self.annot.xy = pos
        idx = np.abs(self.resdfa['zlim']-pos[0]) < 1.e-8
        idx &= np.abs(self.resdfa['nsn']-pos[1]) < 1.e-8
        sel = self.resdfa[idx]
        color = sel['color'].values.item()

        text = "{} \n family: {}".format(
            sel['dbName'].values.item(), sel['family'].values.item())

        self.annot.set_text(text)
        # self.annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
        self.annot.get_bbox_patch().set_facecolor(color)
        self.annot.get_bbox_patch().set_alpha(0.4)

    def hover(self, event):
        """
        Method to plot annotation when pointer is on a point

        Parameters
        ---------------
        event: matplotlib.backend_bases.MouseEvent
          mouse pointer?
        """
        vis = self.annot.get_visible()
        if event.inaxes == self.ax:
            cont, ind = self.sc.contains(event)
            if cont:
                self.update_annot(ind)
                self.annot.set_visible(True)
                self.fig.canvas.draw_idle()
            else:
                if vis:
                    self.annot.set_visible(False)
                    self.fig.canvas.draw_idle()


class plot_DD_Moll:

    def __init__(self, data, dbName, season, nside):

        self.nside = nside
        self.pixel_area = hp.nside2pixarea(nside, degrees=True)

        idx = data['cadence'] == dbName
        idx &= data['season'] == season
        sel = data[idx]

        for fieldName in np.unique(sel['fieldname']):
            ida = sel['fieldname'] == fieldName
            print(fieldName, sel[ida])

        sel = pd.DataFrame(sel)
        sel[['nsn_med_faint', 'zlim_faint']] = sel[[
            'nsn_med_faint', 'zlim_faint']].round(2)
        sel = np.around(sel, decimals=2)
        xmin = 0.01
        xmax = np.max(sel['nsn_med_faint'])
        self.plotMollview(sel, 'nsn_med_faint',
                          'NSN', np.sum, xmin, xmax)

        xmin = 0.01
        xmax = np.max(sel['zlim_faint'])
        self.plotMollview(sel, 'zlim_faint',
                          'zlim', np.median, xmin, xmax)

    def plotMollview(self, data, varName, leg, op, xmin, xmax):
        """
        Method to display results as a Mollweid map

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

        """
        npix = hp.nside2npix(self.nside)

        hpxmap = np.zeros(npix, dtype=np.float)
        hpxmap = np.full(hpxmap.shape, 0.)
        hpxmap[data['healpixID'].astype(
            int)] += data[varName]

        norm = plt.cm.colors.Normalize(xmin, xmax)
        cmap = plt.cm.jet
        cmap.set_under('w')
        resleg = op(data[varName])
        if 'nsn' in varName:
            resleg = int(resleg)
        else:
            resleg = np.round(resleg, 2)
        title = '{}: {}'.format(leg, resleg)
        title = ''
        hp.mollview(hpxmap, min=xmin, max=xmax, cmap=cmap,
                    title=title, nest=True, norm=norm)
        hp.graticule()

        # save plot here
        name = leg.replace(' - ', '_')
        name = name.replace(' ', '_')

        plt.savefig('Plots_pixels/Moll_{}.png'.format(name))


def plotArea(data, nside):
    pixArea = hp.nside2pixarea(nside, degrees=True)

    grpa = area(data, pixArea)
    # select data with zlim_faint>0. and NSN > 10.

    idx = data['zlim_faint'] > 0.
    # idx &= metricValues['nsn_zfaint'] > 10.
    data = data[idx]
    grpb = area(data, pixArea)

    grp = pd.merge(grpa, grpb, on=['cadence',
                                   'dbName_plot'], how='left').reset_index()

    grp['eff_area'] = grp['Npixels_y']/grp['Npixels_x']
    grp = grp.sort_values('cadence')
    figa, axa = plt.subplots()

    axa.plot(grp['dbName_plot'], grp['Npixels_x'],
             color='k', label='total area')
    axa.plot(grp['dbName_plot'], grp['Npixels_y'],
             color='r', label='effective area')
    axa.tick_params(axis='x', labelrotation=20.)
    axa.grid()
    axa.set_ylabel('Median survey area per season [deg$^2$]')
    axa.legend()

    figb, axb = plt.subplots()

    grp = grp.sort_values('eff_area')
    axb.barh(grp['dbName_plot'], grp['eff_area'])
    axb.set_xlabel('Effective area fraction')
    axb.grid()
    # axb.tick_params(axis='x', labelrotation=270.)


def area(data, pixArea):

    grpa = data.groupby(['cadence', 'season', 'fieldname', 'dbName_plot']).apply(
        lambda x: pd.DataFrame({'Npixels': [len(x)*pixArea]})).reset_index()

    grpa_sum = grpa.groupby(['cadence', 'season', 'dbName_plot']).sum()
    grpa_med = grpa_sum.groupby(['cadence', 'dbName_plot']).median()

    return grpa_med
