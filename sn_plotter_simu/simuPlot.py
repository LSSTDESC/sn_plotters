import h5py
import numpy as np
from astropy.table import Table, vstack
import pprint
from . import plt, filtercolors, fontsize
import glob
from sn_tools.sn_utils import multiproc, gather_results
from sn_tools.sn_io import loopStack


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
        parNames = glob.glob(
            '{}/{}/Simu_*{}*.hdf5'.format(self.dbDir, self.dbName, tagName))
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
        fig.suptitle(title, fontsize=thesize)

        for i, var in enumerate(toplot):
            ix = int(i/2)
            iy = i % 2
            axis = ax[ix][iy]
            axis.hist(sel[var], histtype='step')  # bins=len(sel[var]))
            axis.set_xlabel(var, fontsize=20)
            axis.set_ylabel('Number of entries', fontsize=thesize)
            axis.tick_params(axis='x', labelsize=thesize)
            axis.tick_params(axis='y', labelsize=thesize)
        plt.show()

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
