import h5py
# import matplotlib.pyplot as plt
from astropy.table import Table, vstack
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
from . import plt


class FitPlots:
    """
    class to display some results from sn_fit_lc

    Parameters
    ---------------
    files: dict
      data to process
      the dict should have:
        - as keys: labels
        - as values: full path to the files

    """

    def __init__(self, files):

        self.SN_table = {}

        for key, val in files.items():
            self.SN_table[key] = self.load_file(val)

    def load_file(self, fileName):
        """
        Method  to load SN data

        Parameters
        ---------------
        fileName: str
          name of the file to load

        Returns
        -----------
        params: astropy table
        with SN

        """

        f = h5py.File(fileName, 'r')
        print(f.keys(), len(f.keys()))
        sn = Table()
        for i, key in enumerate(f.keys()):
            sni = Table.read(fileName, path=key)
            sn = vstack([sn, sni])

        return sn

    def plot2D(self, tabs, varx, vary, legx, legy):
        """
        2D plot method

        Parameters
        ---------------
        tabs: array
           data to plot
        varx: str
           variable name - x axis
        vary: str
          variable name - y axis
         legx: str
           x-axis legend
        legy: str
          y-axis legend

        """

        fig, ax = plt.subplots()

        self.plot2D_indiv(ax, tabs, varx, vary)

        ax.grid()
        ax.set_xlabel(legx)
        ax.set_ylabel(legy)
        ax.set_ylim([0., 0.1])
        ax.set_xlim([0.01, 1.2])
        ax.legend(loc='upper left')

    def plot2D_indiv(self, ax, tabs, varx, vary, label='', color_cut=0.04):
        """
        2D plot method

        Parameters
        ---------------
        ax: matplotlib axis
           ax for the plot
        tabs: array
           data to plot
        varx: str
           variable name - x axis
        vary: str
          variable name - y axis
        label: str,opt
           label for the display (default: '')
        color_cut: float,opt
           color cut to apply (default: 0.04)

        """

        dict_interp = {}
        for key, tab in tabs.items():

            idx = tab[vary] > 0
            sel = tab[idx]

            ax.plot(sel[varx], np.sqrt(sel[vary]), label=key)

            """
            interp = interp1d(
                np.sqrt(sel[vary]), sel[varx], bounds_error=False, fill_value=0.)

            dict_interp[key] = interp1d(sel[varx], np.sqrt(
                sel[vary]), bounds_error=False, fill_value=0.)

            zlim = interp(color_cut)
            ax.plot(ax.get_xlim(), [color_cut]*2,
                    linestyle='--', color='k')
            ax.plot([zlim]*2, [0., 0.08], linestyle='--', color='k')
            mystr = 'z$_{lim}$'
            ax.text(zlim-0.03, 0.085, '{}={}'.format(mystr, np.round(zlim, 2)))
            """
        # Compare variation
        if len(tabs) < 1:

            zplot = np.arange(0.05, 0.8, 0.01)
            df = pd.DataFrame(zplot.tolist(), columns=['z'])
            for key, val in dict_interp.items():
                kk = '_'.join(key.split('_')[:-1])
                print('kkkk', kk)
                df['sigC_{}'.format(kk)] = val(zplot)

            print(df)
            figb, axb = plt.subplots()
            axb.plot(df['z'], df['sigC_fast_cosmo'] /
                     df['sigC_cosmo_cosmo'], label='fc/cc')
            axb.plot(df['z'], df['sigC_fast_fast'] /
                     df['sigC_cosmo_cosmo'], label='ff/cc')
            axb.plot(df['z'], df['sigC_fast_fast'] /
                     df['sigC_fast_cosmo'], label='ff/fc')

            axb.legend()
            axb.grid()
            axb.set_xlabel('z')
            axb.set_ylabel('$\sigma_{Color}$ ratio')
            axb.set_ylim([0.95, 1.1])
            axb.set_xlim([0.01, 0.78])
