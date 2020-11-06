import numpy as np
import pandas as pd
from sn_tools.sn_io import loopStack
import glob
import healpy as hp

from . import plt
import matplotlib
import matplotlib.patches as mpatches
import tkinter as tk
from tkinter import font as tkFont
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
matplotlib.use('tkagg')


def plot_DDSummary(metricValues, forPlot, sntype='faint'):
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
    sntype: str
      type of the supernova (faint or medium) (default: faint) for display


    Returns
    -----------
    Plot (NSN, zlim)


    """

    # select data with zlim_faint>0. and NSN > 10.

    idx = metricValues['zlim_faint'] > 0.
    # idx &= metricValues['nsn_zfaint'] > 10.
    sel = metricValues[idx]

    print('selec', sel[['zlim_faint', 'nsn_zfaint', 'nsn_med_zfaint']])
    # estimate some stats to display

    data = pd.DataFrame(np.copy(sel))

    summary = data.groupby(['cadence']).agg({'nsn_zfaint': 'sum',
                                             'nsn_med_zfaint': 'sum',
                                             'nsn_zmedium': 'sum',
                                             'zlim_faint': 'median',
                                             'zlim_medium': 'median', }).reset_index()

    summary_fields = data.groupby(['cadence', 'fieldname']).agg({'nsn_zfaint': 'sum',
                                                                 'nsn_med_zfaint': 'sum',
                                                                 'nsn_zmedium': 'sum',
                                                                 'zlim_faint': 'median',
                                                                 'zlim_medium': 'median', }).reset_index()

    summary_fields_seasons = data.groupby(['cadence', 'fieldname', 'season']).agg({'nsn_zfaint': 'sum',
                                                                                   'nsn_med_zfaint': 'sum',
                                                                                   'nsn_zmedium': 'sum',
                                                                                   'zlim_faint': 'median',
                                                                                   'zlim_medium': 'median', }).reset_index()

    # change some of the type for printing
    summary.round({'zlim_faint': 2, 'zlim_medium': 2})
    summary['nsn_zfaint'] = summary['nsn_zfaint'].astype(int)
    summary['nsn_zmedium'] = summary['nsn_zmedium'].astype(int)

    # plot the results

    # per field and per season
    plotNSN(summary_fields_seasons, forPlot, sntype=sntype)
    # per field, for all seasons
    plotNSN(summary_fields, forPlot, sntype=sntype)
    # Summary plot: one (NSN,zlim) per cadence (sum for NSN, median zlim over the fields/seasons)
    plotNSN(summary, forPlot, sntype=sntype)


def plotNSN(summary, forPlot, sntype='faint'):
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
    varx = 'zlim_{}'.format(sntype)
    vary = 'nsn_med_z{}'.format(sntype)
    xshift = 1.0
    yshift = 1.01

    for group in np.unique(forPlot['group']):
        idx = forPlot['group'] == group
        sel = forPlot[idx]
        # print(group, sel['dbName'])
        marker = sel['marker'].unique()[0]
        color = sel['color'].unique()[0]

        print('ici', sel['dbName'].str.strip(), summary['cadence'])
        selcad = summary[summary['cadence'].str.strip().isin(
            sel['dbName'].str.strip())]

        # plot
        ax.plot(selcad[varx], selcad[vary], color=color,
                marker=marker, lineStyle='None')

        # get the centroid of the data and write it
        centroid_x = selcad[varx].mean()
        centroid_y = selcad[vary].mean()
        ax.text(xshift*0.99*centroid_x, yshift *
                1.01*centroid_y, group, color=color, fontsize=fontsize)

    ax.grid()
    ax.set_xlabel('$z_{'+sntype+'}$')
    ax.set_ylabel('$N_{SN} (z<)$')


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
        #self.data = df.loc[:,~df.columns.str.contains('mask', case=False)]
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

        #zlim = self.zlim_med()
        nsn, sig_nsn = self.nSN_tot()
        nsn_extrapol = int(np.round(nsn*self.ratiopixels))

        """
        test = self.data.apply(lambda x: self.ana_filters(x),axis=1,result_type='expand')
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
            resdf['N_{}'.format(ba)] = [med_meds['N_{}'.format(ba)]/resdf['N_total']]
            for bb in 'grizy':
                combi=''.join(sorted('{}{}'.format(ba,bb)))
                #resdf['cadence_{}{}'.format(ba,bb)] = [med_meds['cadence_{}{}'.format(ba,bb)]]
                resdf['N_{}'.format(combi)] = [med_meds['N_{}'.format(combi)]/resdf['N_total']] 
                for bc in 'grizy':
                    combi=''.join(sorted('{}{}{}'.format(ba,bb,bc)))
                    #resdf['cadence_{}{}{}'.format(ba,bb,bc)] = [med_meds['cadence_{}{}{}'.format(ba,bb,bc)]]
                    resdf['N_{}'.format(combi)] = [med_meds['N_{}'.format(combi)]/resdf['N_total']]
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
        #self.ax = ax
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
        self.fig.suptitle('(nSN,zlim) supernovae metric', fontsize=15)
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
        #PlotSummary_Annot(self.resdf, fig=self.fig, ax=self.ax)
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
        #PlotSummary_Annot(resfi, fig=self.fig, ax=self.ax, hlist=highlightlist)
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
        self.fig = plt.Figure(figsize=(14, 8), dpi=100)
        self.ax = self.fig.add_subplot(111)
        leg = 'days$^{-1}$'
        self.fig.suptitle('(nSN,zlim) supernovae metric', fontsize=15)
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
        """
        z_button = tk.Button(button_frame, text="zlim", command=(
            lambda e=ents: self.updateData_nv(e)), bg='yellow', height=heightb, width=widthb, fg='red', font=helv36)
        """
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

        """
        listbox = tk.Listbox(root)
        listbox.insert(0, 'element1')
        listbox.insert(1, 'element2')
        listbox.insert(2, 'element3')
        listbox.pack()
        """
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
        self.plotMetric(resfi, highlightlist)

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

    def plotMetric(self, resdf, hlist=[]):

        self.resdfa = resdf
        x = resdf['zlim'].to_list()
        y = resdf['nsn'].to_list()
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
            y = resh['nsn'].to_list()
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
