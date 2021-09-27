import glob
import h5py
from astropy.table import Table, vstack
from . import plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from sn_tools.sn_rate import SN_Rate


class plotEffi:
    """
    class to plot and process efficiencies from the nsn metric

    Parameters
    --------------
    dbDir: str
      main location dir of the files to process
    dbName: str
      OS to process
    fieldName: str
      DD field to consider

    """

    def __init__(self, dbDir, dbName, fieldName):

        self.data = self.load(dbDir, dbName, fieldName)

    def load(self, dbDir, dbName, fieldName):
        """
        Method to load the data

        Parameters
        --------------
        dbDir: str
        main location dir of the files to process
        dbName: str
          OS to process
        fieldName: str
          DD field to consider

        Returns
        -----------
        astropy table with efficiency data

        """
        fName = '{}/{}/NSN_{}/*.hdf5'.format(dbDir, dbName, fieldName)

        fis = glob.glob(fName)

        fFile = h5py.File(fis[0], 'r')

        keys = list(fFile.keys())
        data = Table()

        for key in keys:
            tab = Table.read(fFile, path=key)
            data = vstack([data, tab])

        print(data)

        return data

    def plotEffi(self):
        """
        Method to plot efficiencies

        """
        for healpixID in np.unique(self.data['healpixID']):
            idx = self.data['healpixID'] == healpixID
            sel = self.data[idx]
            fig, ax = plt.subplots(figsize=(12, 10))
            for season in np.unique(sel['season']):

                self.plotEffi_indiv(sel, ax, healpixID, season, 'effi', 'effi_err',
                                    'Observing Efficiencies', ls='-', label='season {}'.format(season))

            ax.set_xlabel('z')
            ax.set_ylabel('Observing Efficiency')
            ax.legend()
            ax.set_xlim([0.0, 0.7])
            ax.grid()
        plt.show()

    def plotEffi_indiv(self, data, ax, healpixID, season, vary='effi_err', erry=None, legy='', ls='None', sn_x1=-2.0, sn_color=0.2, lineColor='k', label=None, marker='o'):
        """
        Simple method to plot vs z

        Parameters
        --------------
        data: pandas df
          data to plot
        ax:
          axis where to plot
        healpixID: int
          healpixID of the cell to plot
        season: int
          season of the cell to plot
        vary: str
          variable (column of effi) to plot
        erry: str, opt
          error on y-axis (default: None)
        legy: str, opt
          y-axis legend (default: '')
        sn_x1: float, opt
          sn stretch (default: -2.0)
        sn_color: float, opt
          sn color (default: 0.2)
        lineColor: str, opt
         effi curve color (default: black)
        label: str, opt
          curves label (default: None)
        """
        idx = data['healpixID'] == healpixID
        idx &= data['season'] == season
        idx &= data['x1'] == sn_x1
        idx &= data['color'] == sn_color
        grp = data[idx]

        yerr = None

        if erry is not None:
            yerr = grp[erry]

        ax.errorbar(grp['z'], grp[vary], yerr=yerr,
                    marker=marker, label=label, lineStyle=ls, ms=10)

    def getRates(self, survey_area=9.6):

        rateSN = SN_Rate(H0=70., Om0=0.3,
                         min_rf_phase=-15., max_rf_phase=30.)

        # estimate the rates and nsn vs z
        zz, rate, err_rate, nsn, err_nsn = rateSN(zmin=0.01,
                                                  zmax=1.,
                                                  duration=180.,
                                                  survey_area=survey_area,
                                                  account_for_edges=True)

        # rate interpolation
        rateInterp = interp1d(zz, nsn, kind='linear',
                              bounds_error=False, fill_value=0)
        rateInterp_err = interp1d(zz, err_nsn, kind='linear',
                                  bounds_error=False, fill_value=0)

        return rateInterp, rateInterp_err

    def get_zlims(self, data, ax, healpixID, season, sn_x1=-2.0, sn_color=0.2, zlim_coeff=0.95, zplot=np.arange(0.01, 0.7, 0.01)):

        rateInterp, rateInterp_err = self.getRates()

        idx = data['healpixID'] == healpixID
        idx &= data['season'] == season
        idx &= data['x1'] == sn_x1
        idx &= data['color'] == sn_color
        effi = data[idx]

        # interpolate efficiency vs z
        effiInterp = interp1d(
            effi['z'], effi['effi'], kind='linear', bounds_error=False, fill_value=0.)
        # interpolate variance efficiency vs z
        effiInterp_err = interp1d(
            effi['z'], effi['effi_err'], kind='linear', bounds_error=False, fill_value=0.)

        nsn_cum = np.cumsum(effiInterp(zplot)*rateInterp(zplot))
        nsn_cum_err = []
        for i in range(len(zplot)):
            siga = effiInterp_err(zplot[:i+1])*rateInterp(zplot[:i+1])
            sigb = effiInterp(zplot[:i+1])*rateInterp_err(zplot[:i+1])
            nsn_cum_err.append(np.cumsum(
                np.sqrt(np.sum(siga**2 + sigb**2))).item())

        if nsn_cum[-1] >= 1.e-5:
            nsn_cum_norm = nsn_cum/nsn_cum[-1]  # normalize
            nsn_cum_norm_err = nsn_cum_err/nsn_cum[-1]  # normalize
            zlim = interp1d(nsn_cum_norm, zplot,
                            bounds_error=False, fill_value=-1.)
            zlim_plus = interp1d(nsn_cum_norm+nsn_cum_norm_err,
                                 zplot, bounds_error=False, fill_value=-1.)
            zlim_minus = interp1d(
                nsn_cum_norm-nsn_cum_norm_err, zplot, bounds_error=False, fill_value=-1.)
            zlimit = zlim(zlim_coeff).item()
            zlimit_minus = zlim_plus(zlim_coeff).item()
            zlimit_plus = zlim_minus(zlim_coeff).item()

        return nsn_cum_norm, zlimit

    def plotCumul(self, data, ax, healpixID, season, sn_x1=-2.0, sn_color=0.2, zlim_coeff=0.95, label='', ls='None', marker='o'):
        zplot = np.arange(0.01, 0.7, 0.01)
        nsn_cum_norm, zlimit = self.get_zlims(
            data, ax, healpixID, season, sn_x1=sn_x1, sn_color=sn_color, zlim_coeff=zlim_coeff, zplot=zplot)
        ax.plot(zplot, nsn_cum_norm, label=label, ls=ls, marker=marker)
        zcomp = '$z_{complete}^{0.95,'+str(season)+'}$ = '
        ax.text(0.05, 0.958-0.08*season,
                zcomp+'{}'.format(np.round(zlimit, 2)), color=plt.gca().lines[-1].get_color())
        ax.plot([zlimit]*2, [0., 0.95], color='k', ls='dashed')
        # ax.fill_between(zplot, nsn_cum_norm-nsn_cum_norm_err,
        #              nsn_cum_norm+nsn_cum_norm_err, color='y')

    def plotNSN(self, data, ax, healpixID, season, sn_x1=0.0, sn_color=0.0, zlim_coeff=0.95, label='', ls='solid', marker='o'):

        rateInterp, rateInterp_err = self.getRates(survey_area=1.)

        idx = data['healpixID'] == healpixID
        idx &= data['season'] == season
        idx &= data['x1'] == sn_x1
        idx &= data['color'] == sn_color
        effi = data[idx]
        zplot = np.arange(0.01, 0.71, 0.01)
        # interpolate efficiency vs z
        effiInterp = interp1d(
            effi['z'], effi['effi'], kind='linear', bounds_error=False, fill_value=0.)
        # interpolate variance efficiency vs z
        effiInterp_err = interp1d(
            effi['z'], effi['effi_err'], kind='linear', bounds_error=False, fill_value=0.)

        nsn_cum = np.cumsum(effiInterp(zplot)*rateInterp(zplot))

        ax.plot(zplot, nsn_cum, label=label, ls=ls, marker=marker)

        nsn_cum_norm, zlimit = self.get_zlims(
            data, ax, healpixID, season, sn_x1=sn_x1, sn_color=sn_color, zlim_coeff=zlim_coeff, zplot=zplot)
        zplotb = np.arange(0.01, zlimit+0.005, 0.01)
        nsn_zlimit = np.cumsum(effiInterp(zplotb)*rateInterp(zplotb))[-1]
        print(nsn_zlimit)
        #zcomp = '$\mathrm{N_{SN}^{z\leq z_{complete}^{0.95}}}$ = '
        zcomp = '$\mathrm{N_{SN}^{'+str(season)+'}}$='
        if season <= 5:
            xpos = 0.03
            ypos = 14.-2.*season
        else:
            xpos = 0.25
            ypos = 14.-2.*(season-5)

        ax.text(xpos, ypos,
                zcomp+'{}'.format(int(nsn_zlimit)), color=plt.gca().lines[-1].get_color())

        ax.plot([zlimit]*2, [0., nsn_zlimit], ls='dashed', color='k')
        #ax.plot([0., zlimit], [nsn_zlimit]*2, ls='dashed', color='k')
