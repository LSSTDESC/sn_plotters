import pandas as pd
import numpy as np
import numpy.lib.recfunctions as rf
from . import plt


class PlotHist:
    """
    class to plot histograms for a set of OS

    Parameters
    ---------------
    dbDir: str
      loca. directory of the files
    forPlot: pandas df
      configuration parameters

    """

    def __init__(self, dbDir, forPlot):

        self.data = pd.DataFrame()
        for i,row in forPlot.iterrows():
            dbName = row['dbName']
            color = row['color']
            tab = np.load(
                '{}/{}/Global/{}_SNGlobal.npy'.format(dbDir, dbName, dbName))
            df = pd.DataFrame(tab)
            dbName_new = dbName.split('_10yrs')[0]
            df['dbName'] = dbName_new
            df['color'] = color
            df['nights'] = len(tab)
            self.data = pd.concat((self.data, df))

        print('alors',self.data,self.data.columns)

    def plot(self, plotstr, legx, legy='Number of Entries'):
        """
        Method to plot multiple histograms

        Parameters
        ---------------
        plotstr: str
        variable to plot
        legx: str
        x-axis legend
        legy: str, opt
        y-axis legend (default: Number of Entries)

        """

        fig, ax = plt.subplots()
        bins = np.arange(22.,27.,0.1)
        for dbName in np.unique(self.data['dbName']):
            idx = self.data['dbName'] == dbName
            sel = self.data[idx]
            #ib = sel[plotstr]>=0.
            print(dbName,sel[plotstr])
            ax.hist(sel[plotstr], histtype='step', lw=2,
                    label=dbName,density=True)

        ax.legend(ncol=1, loc='best', prop={'size': 12}, frameon=True)
        ax.set_xlabel(legx)
        ax.set_ylabel(legy)


    def plotCumulative(self,plotstr,leg,legx,years=[1]):
        """
        Method to plot the cumulative of plotstr variable for the list of years

        Parameters
        --------------
        plotstr: str
          variable to plot
        years: list(int),opt
          years to consider for plotting (default: 1)

        """

        fig, ax = plt.subplots(figsize=(14,8))
        fig.suptitle(leg)
        for year in years:
            idx = self.data['night'] >= 365.*(year-1)
            idx &= self.data['night'] <365.*year
            sel = self.data[idx]
            for dbName in np.unique(sel['dbName']):
                idxb = sel['dbName'] == dbName
                selb = sel[idxb]
                print(selb.columns)
                t0 = np.min(selb['night'])
                ax.plot(selb['night']-t0,np.cumsum(selb[plotstr]),label=dbName)

        ax.legend()
        ax.set_xlabel('night')
        ax.set_ylabel(legx)

        plt.savefig('Plots/{}_cumul.png'.format(plotstr))
        
    def plotBarh(self, plotstr, legx, legy='Number of Entries'):
        """
        Method to plot multiple histograms

        Parameters
        ---------------
        plotstr: str
        variable to plot
        legx: str
        x-axis legend
        legy: str, opt
        y-axis legend (default: Number of Entries)

        """

        fig, ax = plt.subplots()
        meds = self.data.groupby(['dbName']).median().reset_index()
        ax.barh(meds['dbName'],meds[plotstr])
        """
        for dbName in np.unique(self.data['dbName']):
            idx = self.data['dbName'] == dbName
            sel = self.data[idx]
            print('hhhh',sel.columns)
            ax.barh(sel['dbName'],sel[plotstr])
        """
        ax.legend(ncol=1, loc='best', prop={'size': 12}, frameon=True)
        ax.set_xlabel(legx)
        # ax.set_ylabel(legy)


class PlotTime:
    """
    Class to plot some variable vs time (night number)

    Parameters
    ---------------
    dbDir: str
       location dir of the files
    dbName: str
      OS to display
    forPlot: pandas df
      configuration parameters (marker, color, ...)

    """

    def __init__(self, dbDir, dbName, forPlot):

        self.dbName = dbName

        # loading data
        search_path = '{}/{}/Global/{}_SNGlobal.npy'.format(
            dbDir, dbName, dbName)
        self.data = np.load(search_path)

        # get marker for display
        idx = forPlot['dbName'] == dbName
        self.mark = forPlot[idx]['marker'].values[0]

    def plot(self, varx, legx, vary, legy, nightBeg=0, nightEnd=365):
        """
        Method to plot vary vs varx

        Parameters
        ---------------
        varx: str
          variable to display - xaxis
        vary: str
          variable to display - yaxis
        legx: str
          x-axis legend
        legy: str
          y-axis legend
        nightBeg: int
          first night to consider
        nightEnd: int
          last night to consider


        """

        # select the night
        idx = self.data['night'] < nightEnd
        idx &= self.data['night'] >= nightBeg
        sel = self.data[idx]
        # ax.plot(sel['night'],sel[plotstr],label=dbName,linestyle='',marker='o')

        sel.sort(order='night')

        fig, ax = plt.subplots()

        ax.plot(sel[varx], sel[vary], label=self.dbName,
                ls='None', color='k', marker=self.mark)

        ax.legend(ncol=1, loc='best', prop={'size': 12}, frameon=True)
        ax.set_xlabel(legx)
        ax.set_ylabel(legy)


class PlotStat:
    """
    class to display median values of the Global Metric

    Parameters
    ---------------
    dbDir: str
      location directory of the files
    forPlot: str
      configuration file (dbNames, marker, colors, ...)

    """

    def __init__(self, dbDir, forPlot):

        self.config = forPlot
        # first: estimate stats
        r = []
        res = pd.DataFrame()
        for dbName in forPlot['dbName']:
            tab = np.load(
                '{}/{}/Global/{}_SNGlobal.npy'.format(dbDir, dbName, dbName))
            #rint = [dbName, np.median(tab['nfc']), np.median(tab['obs_area'])]
            df = pd.DataFrame(tab)
            
            df = df.mask(df < 0)
            df['dbName'] = dbName.split('_10yrs')[0]
            df['obs_nights'] = len(tab)
            df['survey_length'] = np.max(tab['night'])-np.min(tab['night'])
            #print(dbName,np.nanmedian(df['med_fiveSigmaDepth_i']))
            meds = df.groupby('dbName').median().reset_index()

            #print('alors',meds['nvisits_u'])
            sums = df.groupby('dbName').sum().reset_index()
            vv = ['dbName']
            for band in 'ugrizy':
                vv.append('frac_{}'.format(band))
                vv.append('nvisits_{}'.format(band))
                sums['frac_{}'.format(band)] = sums['nvisits_{}'.format(band)] /sums['nvisits']
                
            vv.append('nvisits')
            meds = meds.merge(sums[vv], left_on=['dbName'], right_on=['dbName'],suffixes=['','_sum'])
            
            #print('test',meds['med_fiveSigmaDepth_g'])
            res = pd.concat((res,meds))
            
            """    
            dbName_new = dbName.split('_10yrs')[0]
            rint = [dbName_new, np.median(tab['nfc_noddf']),np.median(tab['nfc']),]
            #names = ['dbName', 'nfc_med', 'obs_area_med']
            names = ['dbName', 'nfc_noddf_med','nfc_med']
            for band in 'ugrizy':
                rint += [np.sum(tab['nvisits_{}'.format(band)]) /
                         np.sum(tab['nvisits'])]
                names += ['frac_{}'.format(band)]
            r.append((rint))
            """
        # results stored in a record array
        #self.data = np.rec.fromrecords(r, names=names)
        self.data = res.to_records(index=False)
    def listvar(self):
        """
        Method to list the columns of self.data

        """
        print(self.data.dtype.names)

    def plotBarh(self, plotstr, title,xmin=0.05):
        """
        Method to display results as bar histogram

        Parameters
        ---------------
        plotstr: str
          variable to plot
        title: str
          plot title

        """
        self.data.sort(order=plotstr)
        fig, ax = plt.subplots(figsize=(12,5))
        fig.subplots_adjust(top=0.9,left=0.3)
        #fig, ax = plt.subplots()
        fig.suptitle(title)

        myrange = np.arange(len(self.data))
        ax.barh(myrange, self.data[plotstr])

        for row in self.data:
            print(plotstr,row[['dbName',plotstr]])
        plt.yticks(myrange, self.data['dbName'])
        xmina, xmax = ax.get_xlim()
        xmax = 1.03*np.max(self.data[plotstr])
        ax.set_xlim([xmin, xmax])
        plt.grid(axis='x')
        #plt.tight_layout()
        plt.savefig('Plots/{}.png'.format(plotstr))
