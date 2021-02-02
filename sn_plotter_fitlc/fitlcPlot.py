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

    def plot2D(self, tabs, varx, vary, legx, legy, compare=False):
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

        zlim_max = self.plot2D_indiv(ax, tabs, varx, vary, compare=compare)

        ax.grid()
        ax.set_xlabel(legx)
        ax.set_ylabel(legy)
        ax.set_ylim([0., 0.12])
        ax.set_xlim([0.1, zlim_max+0.1])
        ax.legend(loc='upper left')

    def plot2D_indiv(self, ax, tabs, varx, vary, label='', color_cut=0.04, compare=False):
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

        Returns
        ----------
        zlim_max, float: max zlim value (use for scale display purposes)

        """
        print(tabs.keys())
        dict_interp = {}
        zlims = []
        for key, tab in tabs.items():
            idx = tab[vary] > 0
            idx &= tab['z']>0.1
            sel = tab[idx]
            sel.sort(keys=[varx])
            print(np.unique(sel['z']), sel[varx, vary])
            
            interp = interp1d(
                np.sqrt(sel[vary]), sel[varx], bounds_error=False, fill_value=0.)

            interpv = interp1d(sel[varx], np.sqrt(sel[vary]), bounds_error=False, fill_value=0.)

            dict_interp[key] = interp1d(sel[varx], np.sqrt(
                sel[vary]), bounds_error=False, fill_value=0.)

            #zlim = interp(color_cut)

            zlim = self.zlim(interpv,color_cut)
            zlims.append(zlim)

            ax.plot(sel[varx], np.sqrt(sel[vary]),
                    label='{} - zlim={}'.format(key, np.round(zlim, 2)))

            """
            ax.plot(np.sqrt(sel[vary]),sel[varx], 
                    label='{} - zlim={}'.format(key, np.round(zlim, 2)))
            
            zvals = np.arange(0.2,0.7,0.01)
            colors = np.arange(0.02,0.080,0.001)
            ax.plot(colors,interp(colors),color='k')
            """
            ax.plot(ax.get_xlim(), [color_cut]*2,
                    linestyle='--', color='k')
            """
            ax.plot([zlim]*2, [0., 0.08], linestyle='--', color='k')
            mystr = 'z$_{lim}$'
            """
        # Compare variation

        if compare:

            zplot = np.arange(0.05, 0.8, 0.01)
            df = pd.DataFrame(zplot.tolist(), columns=['z'])
            for key, val in dict_interp.items():
                kk = '_'.join(key.split('_')[:-2])
                print('kkkk', kk)
                kk = key
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

        return np.max(zlims)

    def zlim(self, interp, color_cut):
        """
        Method to estimate zlim

        Parameters
        ---------------
        interp: scipy.interpolate.interp1D
          interpolator to be used (here error_color vs z)
        color_cut: float
          color error used as a reference 

        Returns
        ----------
        z : float
          redshift value corresponding to color_cut

        """

        zvals = np.arange(0.2,1.0,0.005)

        colors = interp(zvals)

       

        ii = np.argmin(np.abs(colors-color_cut))

        print(type(colors),colors[ii],zvals[ii])
        return zvals[ii]
