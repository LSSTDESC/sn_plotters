import h5py
import numpy as np
from astropy.table import Table, vstack
import pprint
from . import plt, filtercolors
import glob
from sn_tools.sn_utils import multiproc, gather_results
from sn_tools.sn_io import loopStack
import pandas as pd
from scipy.interpolate import interp1d
from pandas.plotting import table


class SimuPlot:
    """
    class to analyze and plot simulation output (parameters and LC)

    Parameters
    ---------------
    dbDir: str
      location dir of the files
    dbName: str
       procID of the data to process

    """

    def __init__(self, dbDir, dbName, tagName, nproc=8):

        self.dbDir = dbDir
        self.dbName = dbName
        self.tagName = tagName

        # some display parameters
        self.bands = 'ugrizy'
        self.band_id = dict(
            zip(self.bands, [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]))

        # load simulation parameters
        search_path = '{}/{}/Simu_*{}*.hdf5'.format(
            self.dbDir, self.dbName, tagName)
        print('searching files', search_path)
        parNames = glob.glob(search_path)
        print('Nfiles to load', len(parNames))
        pp = {}
        pp['fichtype'] = 'astropyTable'
        params = multiproc(parNames, pp, self.load_multiproc, nproc)

        #params = self.load_params('{}/{}'.format(self.dbDir, parName))
        # loop on this file using the simuPars list
        # select only LC with status=1
        ik = params['status'] == 1
        self.simuPars = params[ik]

    def load_multiproc(self, data, params={}, j=0, output_q=None):
        """
        Function to load simulation parameters for multiptocessing

        Parameters
        ----------------
        data: list? 
        data to process
        params: dict, opt
        parameters for the multiprocessing

        Returns
        ----------
        loaded data 


        """

        fichtype = params['fichtype']
        dictres = {}
        print('processing', j, len(data))
        for io, fname in enumerate(data):
            res = loopStack([fname], fichtype)
            dictres[io] = res

        finres = gather_results(dictres)

        print('end of proc', j, len(finres))
        if output_q is not None:
            return output_q.put({j: finres})
        else:
            return finres

    def load_params(self, paramFile):
        """
        Function to load simulation parameters

        Parameters
        ---------------
        paramFile: str
          name of the parameter file

        Returns
        -----------
        params: astropy table
        with simulation parameters

        """

        f = h5py.File(paramFile, 'r')
        print(f.keys(), len(f.keys()))
        params = Table()
        for i, key in enumerate(f.keys()):
            pars = Table.read(paramFile, path=key)
            params = vstack([params, pars])

        return params

    def plotParameters(self, season=1, toplot=['x1', 'color', 'daymax', 'z']):
        """ Plot simulation parameters
        parameters ('X1', 'Color', 'DayMax', 'z')

        Parameters
        ---------------
        fieldname: (DD or WFD)
        fieldid: (as given by OpSim)
        tab: recarray of parameters
        season: season

        """

        idx = self.simuPars['season'] == season
        sel = self.simuPars[idx]
        thesize = 15
        ncols = 2
        nrows = int(len(toplot)/2)
        print(nrows, ncols)
        fig, ax = plt.subplots(ncols=ncols, nrows=nrows, figsize=(10, 9))
        #title = '{} - fieldid {} - season {}'.format(fieldname, fieldid, season)
        title = 'season {}'.format(season)
        fig.suptitle(title)

        for i, var in enumerate(toplot):
            ix = int(i/2)
            iy = i % 2
            axis = ax[ix][iy]
            axis.hist(sel[var], histtype='step')  # bins=len(sel[var]))
            axis.set_xlabel(var)
            axis.set_ylabel('Number of entries')
            axis.tick_params(axis='x', labelsize=thesize)
            axis.tick_params(axis='y', labelsize=thesize)
        # plt.show()

    def plotLoopLC(self, pause_time=5):
        """
        Function to plot LC in loop

        Parameters
        ---------------
        pause_time: int, opt
          time of the window persistency (in sec) (default: 5 sec)
        """

        # get all simu files

        simu_path = '{}/{}/Simu_*{}*.hdf5'.format(
            self.dbDir, self.dbName, self.tagName)

        simuFiles = glob.glob(simu_path)

        # loop on simuFiles
        for simuName in simuFiles:
            # get corresponding LC file
            namespl = simuName.split('/')
            namespl[-1] = namespl[-1].replace('Simu', 'LC')
            lcName = '/'.join(namespl)
            print('here', simuName, lcName)
            # loop on simu parameters
            for par in self.load_params(simuName):
                print('status', par['status'])
                lc = Table.read(lcName, path='lc_{}'.format(par['index_hdf5']))
                print('lc', lc.columns)
                self.plotFig(lc, pause_time=pause_time)

        """
        # get LC file
        lcFile = '{}/{}/LC_*{}*.hdf5'.format(self.dbDir,
                                             self.dbName, self.tagName)
        f = h5py.File(lcFile, 'r')
        print(f.keys(), len(f.keys()))

        simpars = self.simuPars

        # for i, key in enumerate(f.keys()):
        for par in simpars:
         """

    def plotFig(self, lc, pause_time):
        """
        Method to plot lc on fig

        Parameters
        ---------------
        lc: astropy table
          lc to plot
        pause_time: float
          time for plot persistency (in secs)

        """
        fig, ax = plt.subplots(ncols=2, nrows=3, figsize=(12, 8))
        pprint.pprint(lc.meta)  # metadata
        figtitle = '($x_1,c$)=({},{})'.format(
            lc.meta['x1'], lc.meta['color'])
        figtitle += ' - z={}'.format(np.round(lc.meta['z'], 2))
        figtitle += ' \n daymax={}'.format(np.round(lc.meta['daymax'], 2))
        fig.suptitle(figtitle)

        # print(lc)  # light curve points
        self.plotLC(lc, ax, self.band_id)
        plt.draw()
        plt.pause(pause_time)
        plt.close()

    def plotLC(self, table, ax, band_id):
        """
        Method to plot produced LC

        Parameters
        ---------------
        table: astropy table
          LC to display
        ax: matplotlib axis
         to display
        band_id: int
         id of the band
        """

        fontsize = 15.
        for band in 'ugrizy':
            i = band_id[band][0]
            j = band_id[band][1]
            # ax[i,j].set_yscale("log")
            idx = table['band'] == 'LSST::'+band
            sel = table[idx]

            ax[i, j].errorbar(sel['time'], sel['flux_e_sec'], yerr=sel['flux_e_sec']/sel['snr_m5'],
                              markersize=200000., color=filtercolors[band], linewidth=1)
            if i > 1:
                ax[i, j].set_xlabel('MJD [day]', {'fontsize': fontsize})
            ax[i, j].set_ylabel('Flux [pe/sec]', {'fontsize': fontsize})
            ax[i, j].text(0.1, 0.9, band, horizontalalignment='center',
                          verticalalignment='center', transform=ax[i, j].transAxes)

    def checkLC(self):
        # get LC file
        lcFile = '{}/LC_{}.hdf5'.format(self.dbDir, self.dbName)
        f = h5py.File(lcFile, 'r')
        print(f.keys(), len(f.keys()))

        # stack all LCs
        lctot = Table()
        ptime = []
        for i, key in enumerate(f.keys()):
            lc = Table.read(lcFile, path=key)
            print(lc.columns)
            ptime.append(lc.meta['ptime'])
            lctot = vstack([lctot, lc])
            # break

        # print(lctot.columns)
        toplot = {}
        if 'gamma_interp' in lctot.columns:
            toplot = dict(zip(['snr_m5', 'gamma', 'flux_e_sec'], [
                'snr_m5_interp', 'gamma_interp', 'flux_e_sec_interp']))

        for key, vv in toplot.items():
            fig, ax = plt.subplots()

            ax.hist(lctot[key]/lctot[vv], histtype='step', bins=100)
            #ax.plot(lctot['band'], lctot[key]/lctot[vv], 'ko')

            ax.set_xlabel('{} ratio'.format(key))
            ax.set_ylabel('number of entries')

        fig, ax = plt.subplots()
        ax.hist(ptime, histtype='step', bins=100)
        print('ptime', np.median(ptime))

        plt.show()

    def plotLoopLC_errmod(self):
        """
        Function to plot LC in loop

        Parameters
        ---------------
        pause_time: int, opt
          time of the window persistency (in sec) (default: 5 sec)
        """

        # get all simu files

        simu_path = '{}/{}/Simu_*{}*.hdf5'.format(
            self.dbDir, self.dbName, self.tagName)

        simuFiles = glob.glob(simu_path)

        # loop on simuFiles
        #fig, ax = plt.subplots(ncols=2, nrows=3, figsize=(12, 8))
        # fig.subplots_adjust(wspace=0.5,hspace=0.5)
        # fig.tight_layout()
        for simuName in simuFiles:
            # get corresponding LC file
            namespl = simuName.split('/')
            namespl[-1] = namespl[-1].replace('Simu', 'LC')
            lcName = '/'.join(namespl)
            print('here', simuName, lcName)
            # loop on simu parameters
            tt = Table()
            for par in self.load_params(simuName):
                #print('status', par['z','status'])
                lc = Table.read(lcName, path='lc_{}'.format(par['index_hdf5']))
                lc['z'] = par['z']
                lc['lambdabar_z'] = lc['lambdabar']/(1.+lc['z'])
                lc['fluxerr_model_rel'] = lc['fluxerr_model']/lc['flux']

                # print('lc',lc.columns)
                idx = lc['snr_m5'] >= 1
                idx &= lc['phase'] >= -10.
                idx &= lc['phase'] <= 50
                tt = vstack([tt, lc[idx]], metadata_conflicts='silent')

            tt.convert_bytestring_to_unicode()
            print('here', tt['band'])
            self.plotFig_cumul(tt.to_pandas(), self.band_id)

    def plotFig_cumul(self, table, band_id):
        """
        Method to plot produced LC

        Parameters
        ---------------
        table: astropy table
          LC to display
        ax: matplotlib axis
         to display
        band_id: int
         id of the band
        """
        table = table.round({'lambdabar': 2, 'lambdabar_z': 2})

        for band in 'grizy':
            fig, ax = plt.subplots(figsize=(8, 8))
            # fig.subplots_adjust(right=0.75)
            self.plot_band(ax, band, table)

    def plotFig_cumul_deprecated(self, ax, table, band_id):
        """
        Method to plot produced LC

        Parameters
        ---------------
        table: astropy table
          LC to display
        ax: matplotlib axis
         to display
        band_id: int
         id of the band
        """
        table = table.round({'lambdabar': 2, 'lambdabar_z': 2})

        for band in 'ugrizy':
            i = band_id[band][0]
            j = band_id[band][1]
            # ax[i,j].set_yscale("log")
            idx = table['band'] == 'LSST::'+band
            sel = table[idx]
            if len(sel) > 0:
                zbins = np.arange(0.01, 1.3, 0.05)
                lambdabar = np.mean(sel['lambdabar'])
                print('ere', band, lambdabar)
                bins = lambdabar/(1.+zbins)
                bins.sort()
                group = sel.groupby(pd.cut(sel.lambdabar_z, bins))
                lambdabar_z = group['lambdabar_z'].median()
                fluxerr = group['fluxerr_model_rel'].median()
                zmeds = group['z'].median()
                ax[i, j].plot(lambdabar_z, fluxerr,
                              color=filtercolors[band], linewidth=1, marker='.')
                #ymin, ymax = ax[i,j].get_ylim()
                # ax[i,j].plot([380.]*2,[ymin,ymax],ls='dashed',color='k')
                # ax[i,j].plot([800.]*2,[ymin,ymax],ls='dashed',color='k')

                ax2 = ax[i, j].twiny()
                ax2.plot(lambdabar_z, fluxerr, linewidth=1,
                         marker='.', color=filtercolors[band])
                df = pd.DataFrame(ax2.get_xticks(), columns=['ticks'])
                df['ticks'] = lambdabar/df['ticks']-1
                df = df.round({'ticks': 1})
                df['ticks'] = df['ticks'].astype(str)
                print('iii', df)
                # ax2.set_xticks(range(len(df['ticks'])),df['ticks'].to_list())
                ax2.set_xticklabels(df['ticks'].to_list())
                ax2.set_xlabel(r'$z$', fontdict={
                               'fontsize': fontsize}, loc='left')

                # ax2.set_ylim((ymin,ymax))

                """
                ymin, ymax = ax[i,j].get_xlim()
                xmin = lambdabar/ymax-1
                xmax = lambdabar/ymin-1
                print('hello',xmin,xmax)
                ax2.set_xlim((xmin,xmax))
                
                ax2.plot([],[])
                """
                #ax2.plot(zmeds,fluxerr,color='k', linewidth=1,marker='.')
                # ax2.invert_xaxis()

                """
                ax2.plot(lambdabar_z,fluxerr,color='k', linewidth=1,marker='.')
                ticks = lambdabar/ax2.get_xticks()-1
                print('ticks',ticks)
                ax2.set_xticks(ticks)
                """
                # ax2.invert_xaxis()
                """
                import matplotlib as mpl
                print(mpl.__version__)
                ax2 = ax[i,j].secondary_xaxis('top', functions=(self.ff,self.inv))
                ax2.invert_xaxis()
                """

            #ax[i, j].set_xlabel(r'$\frac{\bar{\lambda}}{1+z}$ [nm]', {'fontsize': fontsize})
            #ax[i, j].set_xlabel(r'$\bar{\lambda}/(1+z)$ [nm]', fontsize=fontsize,loc='right')
            ax[i, j].set_ylabel(r'flux error model')

            ax[i, j].text(0.1, 0.9, band, horizontalalignment='center',
                          verticalalignment='center', transform=ax[i, j].transAxes)
            ax[i, j].grid()

    def ff(self, x, lambdabar):
        return lambdabar/x-1.

    def inv(x, lambdabar):
        return lambdabar/(1.+x)

    def plot_band(self, ax, band, table):

        fontsize = 15
        idx = table['band'] == 'LSST::'+band
        sel = table[idx]
        z_blue = {}
        lambda_vals = [380, 370, 360, 350]
        r = []
        if len(sel) > 0:
            zbins = np.arange(0.01, 1.3, 0.05)
            lambdabar = np.mean(sel['lambdabar'])
            for vv in lambda_vals:
                z_blue[vv] = np.round(lambdabar/vv-1, 2)
                r.append((vv, np.round(lambdabar/vv-1, 2)))
            print('hello', band, z_blue)
            bins = lambdabar/(1.+zbins)
            bins.sort()
            group = sel.groupby(pd.cut(sel.lambdabar_z, bins))
            lambdabar_z = group['lambdabar_z'].median()
            fluxerr = group['fluxerr_model_rel'].mean()
            fluxerr_std = group['fluxerr_model_rel'].std()
            zmeds = group['z'].median()
            interp = interp1d(
                zmeds, fluxerr, bounds_error=False, fill_value=0.)
            interpp = interp1d(zmeds, fluxerr+fluxerr_std,
                               bounds_error=False, fill_value=0.)
            interpm = interp1d(zmeds, fluxerr-fluxerr_std,
                               bounds_error=False, fill_value=0.)
            df = pd.DataFrame(r, columns=['blue_cutoff', 'z_blue'])
            df['error'] = np.round(interp(df['z_blue']), 3)
            df['error_plus'] = np.round(interpp(df['z_blue']), 3)
            df['error_minus'] = np.round(interpm(df['z_blue']), 3)
            print(df)
            #table(ax, df, loc="upper right", colWidths=[0.2]*5);
            columns = df.columns.to_list()
            columns = ['blue cutoff [nm]',
                       '$z_{cutoff}$', '$\sigma$', '$\sigma+$', '$\sigma-$']
            #rows = df['blue_cutoff'].to_list()
            cell_text = []
            for j, vals in df.iterrows():
                cell_text.append((vals['blue_cutoff'], vals['z_blue'],
                                 vals['error'], vals['error_plus'], vals['error_minus']))

            the_table = plt.table(cellText=cell_text,
                                  # rowLabels=rows,
                                  colLabels=columns,
                                  loc='upper right',
                                  colWidths=[0.18]*5,
                                  cellLoc='center')
            # bbox=(1., 0.5,0.5,0.5))

            #ax.plot(lambdabar_z,fluxerr,color=filtercolors[band], linewidth=1,marker='.')
            ax.errorbar(lambdabar_z, fluxerr, yerr=fluxerr_std,
                        color=filtercolors[band], linewidth=1, marker='.')
            ax2 = ax.twiny()
            # ax2.plot(lambdabar_z,fluxerr,linewidth=1,marker='.',color=filtercolors[band])
            ax2.errorbar(lambdabar_z, fluxerr, yerr=fluxerr_std,
                         color=filtercolors[band], linewidth=1, marker='.')
            df = pd.DataFrame(ax2.get_xticks(), columns=['ticks'])
            df['ticks'] = lambdabar/df['ticks']-1
            df = df.round({'ticks': 1})
            df['ticks'] = df['ticks'].astype(str)
            ax2.set_xticklabels(df['ticks'].to_list())
            ax2.set_xlabel(r'$z$')

            ax.set_xlabel(r'$\frac{\bar{\lambda}}{1+z}$ [nm]', fontsize=20)
            ax.set_ylabel(
                r'$\sigma_{\mathrm{flux}}^{\mathrm{error \, model}}$/flux')

            ax.text(0.9, 0.1, '{}-band'.format(band), horizontalalignment='center',
                    verticalalignment='center', transform=ax.transAxes)
            ax.grid()
