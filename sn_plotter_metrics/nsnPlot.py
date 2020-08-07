import numpy as np
import pandas as pd
from sn_tools.sn_io import loopStack
import glob
import healpy as hp

from . import plt


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


class NSNAnalysis:
    def __init__(self, dbDir, dbInfo,
                 metricName='NSN', fieldType='WFD',
                 nside=64,
                 x1=-2.0, color=0.2, npixels=-1):
        """
        class to analyze results from NSN metric

        Parameters
        ---------------
        dbDir: str
          location directory where the files to process are
        dbInfo: pandas df
          info from observing strategy (dbName, plotName, color, marker)
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

        self.sntype = 'faint'
        if x1==0. and color==0.:
            self.sntype='medium'
        
        self.dbInfo = dbInfo

        # loading data (metric values)
        search_path = '{}/{}/{}/*NSNMetric_{}*_nside_{}_*.hdf5'.format(
            dbDir, dbInfo['dbName'], metricName, fieldType, nside)
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
        self.data = self.data.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)
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
        bandstat = ['u','g','r','i','z','y','gr','gi','gz','iz','uu','gg','rr','ii','zz','yy']
        
        for b in bandstat:
            self.data['cadence_{}'.format(
                b)] = self.data['season_length']/self.data['N_{}'.format(b)]
            
        meds = self.data.groupby(['healpixID']).median().reset_index()
        meds = meds.round({'zlim_{}'.format(self.sntype): 2})
        med_meds = meds.median()

        resdf = pd.DataFrame([med_meds['zlim_{}'.format(self.sntype)]], columns=['zlim'])
        resdf['nsn'] = [nsn]
        resdf['sig_nsn'] = [sig_nsn]
        resdf['nsn_extra'] = [nsn_extrapol]
        resdf['dbName'] = self.dbInfo['dbName']
        resdf['plotName'] = self.dbInfo['plotName']
        resdf['color'] = self.dbInfo['color']
        resdf['marker'] = self.dbInfo['marker']
        resdf['cadence'] = [med_meds['cadence']]
        resdf['season_length'] = [med_meds['season_length']]
        resdf['N_total'] = [med_meds['N_total']]
        
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

    def ana_filters(self,grp):
        
        print('io',grp['N_filters_night'])
        filter_night = grp['N_filters_night'].split('/')
        print(filter_night)
        r = []
        for vv in filter_night:
            if vv != '':
                vs = vv.split('*')
                print('what',vs)
                r.append((int(vs[0]),vs[1]))

        print(r)
        res = pd.DataFrame(r, columns=['Nvisits','filter'])

        print(res['Nvisits'].sum(),grp['N_total'])
        
        
        return pd.DataFrame({'test':grp})
    
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
        name = leg.replace(' - ','_')
        name = name.replace(' ','_')

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
