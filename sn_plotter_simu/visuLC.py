from sn_tools.sn_io import Read_LightCurve, get_meta, load_SN
from sn_fitter.fit_sn_cosmo import Fit_LC
from astropy.table import Table, vstack


class VisuLC:
    def __init__(self, metaDir, metaFile,
                 SNFile=None, SNDir=None):
        """
        Class to visualize (and fit) LCs

        Parameters
        ----------
        metaFileInput : str
            metadata file.
        metaDirInput : str
            location dir of meta data file.
        SNFileInput : str, optional
             SN file. The default is None.
        SNDirInput : TYPE, optional
             DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """

        meta = Read_LightCurve(file_name=metaFile, inputDir=metaDir)

        paths = meta.get_path()

        self.lcs = {}
        self.metaTot = meta.get_all_data()

        """
        meta = Read_LightCurve(file_name=metaFileInput, inputDir=metaDirInput)

        paths = meta.get_path()

        self.lcs = {}
        self.metaTot = meta.get_all_data()
        """
        """
        for pp in paths:
            if 'table_column_meta' in pp:
                continue
            metaTable = meta.get_table(path=pp)

            metadata = metaTable.meta

            # get lc
            lcDir = metadata['lc_dir']
            lcName = metadata['lc_fileName']
            print('check here', lcName, lcName.replace('LC', 'SN'))

            self.lcs[pp] = Read_LightCurve(file_name=lcName, inputDir=lcDir)

            # print SNIDS
            metaTable['path'] = pp
            metaTot = vstack([metaTot, metaTable])

        print(metaTot['SNID'])
        self.metaTot = metaTot
        """
        print(self.metaTot['SNID'])
        # fit instance
        self.fit = Fit_LC(model='salt3', version='2.0')

        # getting SN (if any)
        self.SN = Table()
        if SNFile != 'None':
            from sn_tools.sn_io import loopStack
            path = '{}/{}'.format(SNDir, SNFile)
            self.SN = loopStack([path], 'astropyTable')

            # print(self.SN.columns, len(self.SN))

    def plot(self, lcpath):
        """
        Method to plot+fit LC corresponding to lcpath

        Parameters
        ----------
        lcpath : str
            path to the LC.

        Returns
        -------
        None.

        """

        idx = self.metaTot['SNID'] == lcpath

        metadata = self.metaTot[idx]
        print(metadata)
        # get lc
        lcDir = metadata['lc_dir'].value[0]
        lcName = metadata['lc_fileName'].value[0]

        lcs = Read_LightCurve(file_name=lcName, inputDir=lcDir)

        lc = lcs.get_table(lcpath)

        # trying to fit here
        outfit = self.fit(lc)

        if len(self.SN) > 0:
            self.plot_SN(lcpath, lc)

    def plot_SN(self, lcpath, lc):
        """
        Method to plot SN info and LC tagged by lcpath

        Parameters
        ----------
        lcpath : str
            SN id.
        lc : atropy table
            corresponding light curve.

        Returns
        -------
        None.

        """

        idl = self.SN['SNID'] == lcpath
        SNsel = self.SN[idl]
        ll = ['SNID', 'x1', 'color', 'daymax', 'n_epochs_m10_p35',
              'n_epochs_m10_p5', 'n_epochs_p5_p20', 'n_bands_m8_p10',
              'selected']

        SNsel.round({'x1': 2, 'color': 4, 'daymax': 1})
        print(SNsel[ll])
        import matplotlib.pyplot as plt
        idx = lc['flux']/lc['fluxerr'] >= 1
        sel_lc = lc[idx]

        fig, ax = plt.subplots(figsize=(7, 9))
        import numpy as np
        colors = dict(zip('ugrizy', ['b', 'c', 'g', 'y', 'r', 'm']))
        for band in np.unique(sel_lc['band']):
            color = colors[band[-1]]
            ido = sel_lc['band'] == band
            sel_b = sel_lc[ido]
            ax.errorbar(sel_b['phase'], sel_b['flux'],
                        yerr=sel_b['fluxerr'],
                        marker='o', color=color, ls='None', label=band)
        ax.grid()
        ax.set_ylabel('flux (pe/s)')
        ax.set_xlabel('phase')
        plt.legend()
        plt.show(block=False)


class SNToLC:
    def __init__(self, metaDir,
                 SNFile, SNDir):
        """
        class to link SN to its LC

        Parameters
        ----------
        metaDir : str
            metadata (simu) dir.
        SNFile : str
            file for SN.
        SNDir : str
            dir for SN.

        Returns
        -------
        None.

        """

        # fit instance
        self.fit = Fit_LC(model='salt3', version='2.0', outType='dict_res')

        # load SN
        SN = load_SN(SNDir, SNFile)

        # get production ID
        prodID = SNFile.split('.hdf5')[0]
        prodID = prodID[3:]

        # get corresponding simu meta data
        meta = get_meta(prodID, metaDir)

        self.sn_vs_lc(SN, meta)

    def sn_vs_lc(self, SN, meta):
        """
        Method to display LC corresponding to SN 

        Parameters
        ----------
        SN : astropy table
            SN.
        meta : astropy table
            metadata.

        Returns
        -------
        None.

        """

        lcs = {}
        io = 0

        outdir = 'OutFig'
        from sn_tools.sn_io import checkDir
        checkDir(outdir)
        for vv in SN:
            io += 1
            # if io >= 3:
            #    continue

            snid = vv['SNID']
            idx = meta['SNID'] == snid
            metadata = meta[idx]
            print(metadata)
            # get lc
            lcDir = metadata['lc_dir'].value[0]
            lcName = metadata['lc_fileName'].value[0]
            if not lcs or lcName not in lcs.keys():
                lcs[lcName] = Read_LightCurve(file_name=lcName, inputDir=lcDir)

            lc = lcs[lcName].get_table(snid)

            # trying to fit here
            outfit = self.fit(lc)
            import matplotlib.pyplot as plt
            fig = self.fit.plotIt(outfit['lc'],
                                  outfit['fitted_model'],
                                  outfit['res_errors'],
                                  outfit['fitstatus'])
            outname = '{}/SN_{}.png'.format(outdir, io)
            fig.savefig(outname)
            plt.close(fig)
