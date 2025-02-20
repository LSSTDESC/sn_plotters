from sn_tools.sn_io import Read_LightCurve, get_meta, load_SN
# from sn_fitter.fit_sn_cosmo import Fit_LC
from astropy.table import Table, vstack
from sn_telmodel.sn_throughputs import get_telescope
import sncosmo
import matplotlib.pyplot as plt
import numpy as np
from random import gauss
import glob
import numpy.lib.recfunctions as rf


class VisuLC:
    def __init__(self, metaDir, metaFile, prodID='None',
                 SNFile='None', SNDir='None',
                 airmassType='const', tag_tel='1.9',
                 airmass=1.2,
                 aerosol=0.0,
                 pwv=4.0,
                 oz=400,
                 remove_sat=0):
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
        SNDirInput : str, optional
             SN input dir. The default is None.
        airmassType: str, optional.
             airmass type for LC fit. The default is const.
        tag_tel: str, opt
          tag for telescope version. The default is 1.9
        remove_sat: int, opt
          To remove LC saturated points when fitting.

        Returns
        -------
        None.

        """

        meta = get_meta(prodID, metaFile, metaDir)

        print('meta', meta)

        """
        paths = meta.get_path()

        self.lcs = {}
        self.metaTot = meta.get_all_data()
        """

        self.metaTot = meta
        print('passed here')

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

        self.metaTot['z'] = np.round(self.metaTot['z'], 4)
        self.metaTot['SNID', 'z'].pprint_all()

        # grab a telescope
        tel_dir = 'throughputs'
        throughputsDir = 'baseline'
        atmosDir = 'atmos'

        telb = '{}_{}'.format(tel_dir, tag_tel)
        through_dir = '{}/{}'.format(telb, throughputsDir)
        atmos_dir = '{}/{}'.format(telb, atmosDir)
        telescope = get_telescope(tel_dir=telb,
                                  through_dir=through_dir,
                                  atmos_dir=atmos_dir,
                                  tag=tag_tel, airmass=airmass,
                                  aerosol=aerosol, pwv=pwv, oz=oz)

        # fit instance
        # self.fit = Fit_LC(model='salt3', version='2.0', telescope=telescope)
        self.prepare_fit(telescope, model='salt3', airmassType=airmassType,
                         airmass=airmass,
                         aerosol=aerosol,
                         pwv=pwv,
                         oz=oz)

        # getting SN (if any)
        self.SN = Table()
        if SNFile != 'None':
            from sn_tools.sn_io import loopStack
            path = '{}/{}'.format(SNDir, SNFile)
            self.SN = loopStack([path], 'astropyTable')

        self.remove_sat = remove_sat
        # print(self.SN.columns, len(self.SN))

    def prepare_fit(self, telescope,
                    model='salt2-extended', version='2.0',
                    airmassType='const', airmass=1.2,
                    aerosol=0.0,
                    pwv=4.0,
                    oz=400):
        """
        Method to load tel bandpasses for sncosmo

        Parameters
        ----------
        telescope : sn_telescope
            Telescope.
        model : str, optional
            Fitter model. The default is 'salt3'.
        version : str, optional
            Fitter version. The default is '2.0'.
        airmassType: str, optional.
            airmass type. The default is 'const'

        Returns
        -------
        None.

        """

        self.register_bands(telescope, airmass=airmass,
                            aerosol=aerosol,
                            pwv=pwv,
                            oz=oz)

        source = sncosmo.get_source(model, version)

        """
        if model == 'salt3':
            source._wave[0] = 1700.
            source._wave[-1] = 24990.
        """
        print('model version', model, version)
        # get the dust
        dustmap = sncosmo.OD94Dust()
        self.model = sncosmo.Model(source=source,
                                   effects=[dustmap, dustmap],
                                   effect_names=['host', 'mw'],
                                   effect_frames=['rest', 'obs'])

    def register_bands(self, telescope, airmass=1.2,
                       aerosol=0.0,
                       pwv=4.0,
                       oz=400):
        """
        Method to register bands in sncosmo

        Returns
        -------
        None.

        """
        from sn_tools.sn_utils import register_bands_sncosmo
        import pandas as pd
        airmass = [airmass]
        pwvs = [pwv]
        ozs = [oz]
        aerosols = [aerosol]

        # values in pandas df
        cols = ['airmass', 'pwv', 'ozone', 'aerosol']
        vals = [airmass, pwvs, ozs, aerosols]
        df = pd.DataFrame.from_dict(dict(zip(cols, vals)))

        """
        print(df_dict)

        df = df_dict['airmass']
        for col in cols[1:]:
            df = df.merge(df_dict[col], how='cross')
        """
        for i, row in df.iterrows():
            airmass = row['airmass']
            aerosol = row['aerosol']
            pwv = row['pwv']
            ozone = row['ozone']
            register_bands_sncosmo(sncosmo,
                                   telescope,
                                   airmass, aerosol, pwv, ozone)

    def plot(self, lcpath):

        if lcpath == 'NSN':
            print('NSN', len(self.metaTot))
        else:
            self.plot_lc(lcpath)

    def plot_lc(self, lcpath):
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
        # print(metadata)
        # get lc
        lcDir = metadata['lc_dir'].value[0]
        lcName = metadata['lc_fileName'].value[0]

        lcs = Read_LightCurve(file_name=lcName, inputDir=lcDir)

        lc = lcs.get_table(lcpath)

        # print('stretch and color', lc.meta['x1'], lc.meta['color'])
        idx = lc['fluxerr'] > 0.
        idx &= lc['flux'] >= 0.
        idx &= lc['snr_m5'] >= 1.
        if self.remove_sat:
            idx &= lc['sat'] == 0
        lc = lc[idx]

        if 'sat' in lc.columns:
            print(lc[['band', 'night', 'flux', 'zp',
                      'zpsys', 'fluxerr', 'seeingFwhmEff', 'sat']])

        # trying to fit here
        # outfit = self.fit(lc)

        # print('fit', outfit)

        # printInfo = len(self.SN) > 0
        # self.plot_SN(lcpath, lc, printInfo)

        sigma_z = 1.e-5
        z = lc.meta['z']
        zmeas = z+gauss(0, sigma_z*(1+z))

        del lc['filter']
        self.model.set(z=lc.meta['z'])

        result, fitted_model = self.fitIt(lc)
        # print(result)
        # print(fitted_model)
        # fitted_model = None

        """
        lc = lc.to_pandas()
        lc['band'] = lc['band'].str.split('_').str.get(0)

        lc = Table.from_pandas(lc)
        print('kkk', lc['z'])
        """
        if fitted_model is not None:
            sncosmo.plot_lc(lc,
                            model=fitted_model,
                            errors=result.errors,
                            zp=25.,
                            xfigsize=8, pulls=False, figtextsize=1.5)
        else:
            sncosmo.plot_lc(lc, xfigsize=9)

        plt.show(block=False)

    def fitIt(self, lc):
        """
        Method to fit a light curve

        Parameters
        ----------
        lc : astropy table
            LC to fit.

        Returns
        -------
        result : array
            fit results.
        fitted_model : array
            fitted model.

        """

        result = None
        fitted_model = None
        try:
            bounds = {'x1': (-3.0, 3.0), 'c': (-0.3, 0.3)}
            self.model.set(mwebv=lc.meta['ebvofMW'])
            # self.model.set(mwebv=0.)
            # print(self.model)
            result, fitted_model = sncosmo.fit_lc(lc,
                                                  self.model,
                                                  vparam_names=[
                                                      't0', 'x0', 'x1', 'c'],
                                                  bounds=bounds,
                                                  minsnr=1.)
        except (RuntimeError, TypeError, NameError) as err:
            print('fit crashed')

        return result, fitted_model

    def plot_SN(self, lcpath, lc, printInfo=False):
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

        if printInfo:
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

        ttext = self.get_text(lc.meta)
        ax.text(0.1, 1.07, ttext, horizontalalignment='center',
                verticalalignment='center',
                transform=ax.transAxes, fontsize=15)

        ttextb = self.get_text(lc.meta,
                               ddict=dict(zip(['z', 'daymax'], [2, 1])))
        ax.text(0.4, 1.07, ttextb, horizontalalignment='center',
                verticalalignment='center',
                transform=ax.transAxes, fontsize=15)

        plt.show(block=False)

    def get_text(self, meta, ddict=dict(zip(['x1', 'color'], [1, 2]))):
        """
        Method to write a text from metadata

        Parameters
        ----------
        meta : dict
            metadata.
        ddict : dict, optional
            what to write (var, round). 
            The default is dict(zip(['x1', 'color'], [1, 2])).

        Returns
        -------
        ttext : str
            the text.

        """

        ttext = ''
        for key, vals in ddict.items():
            ttext += self.simple_text(key, meta, vals)
            ttext += '\n'

        return ttext

    def simple_text(self, var, meta, rounding):
        """
        Method to write a simple text

        Parameters
        ----------
        var : str
            var to write.
        meta : dict
            metadata.
        rounding : int
            rounding for the writing.

        Returns
        -------
        ttext : str
            the text.

        """

        import numpy as np
        ttext = '{}: {}'.format(var, np.round(meta[var], rounding))

        return ttext


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
        self.fit = Fit_LC(model='salt2-extended',
                          version='2.0', outType='dict_res')

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

        idx = SN['n_epochs_bef'] >= 4
        idx &= SN['n_epochs_bef'] >= 10
        SN = SN[idx]

        for vv in SN:
            io += 1
            # if io >= 3:
            #    continue

            if not vv['selected']:
                continue
            snid = vv['SNID']
            idx = meta['SNID'] == snid
            metadata = meta[idx]
            # print(metadata.keys())
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


class VisuNight:
    def __init__(self, dbDir, dbName, fields, colors, colName):
        """
        class to display filter alloc a given night

        Parameters
        ----------
        dbDir : str
            Data location dir.
        dbName : str
            OS to process.
        fields : str
            List of fields to process.
        colors : str
            colors corresponding to fields.
        colName : str
            colName to tag fields.

        Returns
        -------
        None.

        """

        fName = '{}/{}.npy'.format(dbDir, dbName)

        data = np.load(fName, allow_pickle=True)
        # select DDFs
        idx = np.in1d(data[colName], fields.split(','))
        self.ddf = data[idx]

        self.ddplot = dict(zip(fields.split(','), colors.split(',')))
        self.colName = colName

    def plot(self, night):
        """
        Method to plot a night

        Parameters
        ----------
        night : int
            night number.

        Returns
        -------
        None.

        """

        idx = self.ddf['night'] == night
        sel = self.ddf[idx]

        fig, ax = plt.subplots(figsize=(10, 8))
        fig.suptitle('night {}'.format(night))
        fig.subplots_adjust(right=0.80)

        mjd0 = sel['mjd'].min()
        ttime = (sel['mjd']-mjd0)*24.
        sel = rf.append_fields(sel, 'time', ttime)
        for key, vals in self.ddplot.items():
            idxb = sel[self.colName] == key
            selb = sel[idxb]
            if len(selb) > 0:
                ax.plot(selb['time'], selb['filter'], color=vals,
                        linestyle='None', marker='o', label=key.split('DD:')[-1])

        ax.set_xlabel(r'obs time [h]')
        ax.grid(visible=True)
        # ax.legend(bbox_to_anchor=(1., 0.5),
        #          ncol=1, fontsize=12, frameon=False)
        ax.legend(loc='upper left', bbox_to_anchor=(1.0, 0.5),
                  ncol=1, frameon=False, fontsize=15)
        plt.show(block=False)
