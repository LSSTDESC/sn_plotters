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
from sn_tools.sn_utils import register_bands_sncosmo
import pandas as pd
import os

class VisuLC:
    def __init__(self, metaDir, metaFile, snFile='None',prodID='None',
                 airmassType='const', tag_tel='1.9',
                 airmass=1.2,
                 aerosol=0.0,
                 pwv=4.0,
                 ozone=400,
                 remove_sat=0,
                 fit_coadd=0):
        """
        Class to visualize (and fit) LCs

        Parameters
        ----------
        metaFile : str
            metadata file.
        metaDir : str
            location dir of meta data file.
        snFile: str, optional
            SN location dir. The default is None.
        prodID: str, optional
            production ID. The default is None.
        airmassType: str, optional.
             airmass type for LC fit. The default is const.
        tag_tel: str, opt
          tag for telescope version. The default is 1.9
        remove_sat: int, opt
          To remove LC saturated points when fitting.
        fit_coadd: int, opt.
          To fit coadded (band/night) LC points.

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

        #grab the SN
        self.sn_data = pd.DataFrame()
        if snFile != 'None':
            self.sn_data = pd.read_hdf(snFile)
            

        self.metaTot['z'] = np.round(self.metaTot['z'], 4)
        self.metaTot['SNID', 'z'].pprint_all()

        # grab a telescope
        tel_dir = 'throughputs'
        throughputsDir = 'baseline'
        atmosDir = 'atmos'

        telb = '{}_{}'.format(tel_dir, tag_tel)
        through_dir = '{}/{}'.format(telb, throughputsDir)
        atmos_dir = '{}/{}'.format(telb, atmosDir)
        self.telescope = get_telescope(tel_dir=telb,
                                       through_dir=through_dir,
                                       atmos_dir=atmos_dir,
                                       tag=tag_tel, airmass=airmass,
                                       aerosol=aerosol, pwv=pwv, ozone=ozone)

        # fit instance
        # self.fit = Fit_LC(model='salt3', version='2.0', telescope=telescope)
        self.prepare_fit(model='salt3', airmassType=airmassType,
                         airmass=airmass,
                         aerosol=aerosol,
                         pwv=pwv,
                         ozone=ozone)

        self.remove_sat = remove_sat
        # print(self.SN.columns, len(self.SN))

        self.fit_coadd = fit_coadd

    def prepare_fit(self,
                    model='salt2-extended', version='2.0',
                    airmassType='const', airmass=1.2,
                    aerosol=0.0,
                    pwv=4.0,
                    ozone=400):
        """
        Method to load tel bandpasses for sncosmo

        Parameters
        ----------
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

        """
        self.register_bands(telescope, airmass=airmass,
                            aerosol=aerosol,
                            pwv=pwv,
                            oz=oz)
        """
        source = sncosmo.get_source(model, version)

        if model == 'salt3':
            source._wave[0] = 1700.
            source._wave[-1] = 24990.

        print('model version', model, version)
        # get the dust
        dustmap = sncosmo.OD94Dust()
        self.model = sncosmo.Model(source=source,
                                   effects=[dustmap, dustmap],
                                   effect_names=['host', 'mw'],
                                   effect_frames=['rest', 'obs'])

    def register_bands_on_the_fly(self, data):
        """
        Method to register bands on sncosmo

        Parameters
        ----------
        data: pandas df
            data to register

        Returns
        -------
        None.

        """

        print('band registry')
        for i, row in data.iterrows():
            bandname = row['band_cosmo']
            band = row['filter']
            airmass = row['airmass']
            pwv = row['pwv']
            ozone = row['ozone']
            aerosol = row['aerosol']
            register_bands_sncosmo(sncosmo, self.telescope,
                                   bandname, band,
                                   airmass, pwv, ozone, aerosol)

    def register_bands_deprecated(self, telescope, airmass=1.2,
                       aerosol=0.0,
                       pwv=4.0,
                       ozone=400):
        """
        Method to register bands in sncosmo

        Returns
        -------
        None.

        """

        import pandas as pd
        airmass = [airmass]
        pwvs = [pwv]
        ozs = [ozone]
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
        lcDir = metadata['lc_dir'].value[0]
        lcName = metadata['lc_fileName'].value[0]
        
        lc_plus_sn = lc_sn(lcDir,lcName,self.sn_data)
        
        lc_plus_sn.get_infos(lcpath)
        
        lc_plus_sn.plot_all()
        
        return
        print(test)
        
        #
        
        ccols = ['x1','color','daymax','z','x0','ebvofMW']
        ccols_fit = ['x1_fit','color_fit','t0_fit','z_fit','x0_fit','ebvofMW']
        corresp = dict(zip(ccols_fit,ccols))

        pp = {}
        pp_fit={}
        if len(self.sn_data) > 0:
            idx = self.sn_data['SNID'] == lcpath
            sn_fit = self.sn_data[idx]
            print(sn_fit.columns)
            ppa = sn_fit[ccols_fit]
            pp_fit = ppa.to_dict(orient='list')
            ppa = ppa.rename(columns=corresp)
            ppf = ppa.to_dict(orient='list')
            sn_flux_fit = self.grab_fluxes(**ppf)
            pp = sn_fit[ccols].to_dict(orient='list')
            sn_flux_orig = self.grab_fluxes(**pp)
            
        
        # print(metadata)
        # get lc
        lcDir = metadata['lc_dir'].value[0]
        lcName = metadata['lc_fileName'].value[0]

        lcs = Read_LightCurve(file_name=lcName, inputDir=lcDir)

        lc = lcs.get_table(lcpath)

        # get fitted flux here
        sn_fluxes_fit = sn_flux_fit.get_flux(lc)
        sn_fluxes_orig= sn_flux_orig.get_flux(lc)

        lc["phase"] = (lc['time']-pp_fit['t0_fit'][0])/(1.+pp_fit['z_fit'][0])
        filters = np.unique(lc['filter'])
        timescale = 'phase'
    
        for filt in filters:
            idx = lc['filter']==filt
            sel_lc = lc[idx]
            idxb = sn_fluxes_fit['filter'] == filt
            sel_flux_fit = sn_fluxes_fit[idxb]
            idxc = sn_fluxes_orig['filter'] == filt
            sel_flux_orig = sn_fluxes_orig[idxb]
            
            fig, ax = plt.subplots()
            
            ax.errorbar(sel_lc[timescale],sel_lc['flux'],yerr=sel_lc['fluxerr'])
            
            ax.plot(sel_flux_fit[timescale],sel_flux_fit['flux'])
            
            ax.plot(sel_flux_orig[timescale],sel_flux_orig['flux'])
            
            plt.show()
            
            
        
        
        print(test)
        # coadd LC points (band/night) if necessary

        if self.remove_sat:
            idx = lc['sat'] == 0
            lc = lc[idx]

        ccols = ['night', 'airmass', 'ozone', 'aerosol', 'mean_wave', 'band',
                 'pwv', 'zp', 'time', 'band_cosmo', 'zpsys', 'flux', 'fluxerr',
                 'snr_m5', 'snr', 'filter']
        print('before')
        print(lc[['flux', 'fluxerr', 'zp', 'snr_m5', 'snr']])

        if self.fit_coadd:
            mymeta = lc.meta
            df = lc[ccols].to_pandas()
            lc = df.groupby(['filter', 'night']).apply(
                lambda x: self.coadd_lc(x)).reset_index()
            lc['band_cosmo'] = self.telescope.site_name+'::' + \
                lc['filter']+'_' + \
                lc['airmass'].astype(str)+'_' + \
                lc['pwv'].astype(str)+'_' + \
                lc['ozone'].astype(str)+'_' +\
                lc['aerosol'].astype(str)
            lc['band'] = lc['band_cosmo']
            # lc['snr_m5'] = lc['flux']/lc['fluxerr']
            # lc['snr'] = lc['flux']/lc['fluxerr']
            lc = Table.from_pandas(lc)
            lc.meta = mymeta

        print('after')
        print(lc['flux', 'fluxerr', 'zp', 'snr_m5', 'snr'])
        
        # print('stretch and color', lc.meta['x1'], lc.meta['color'])
        idx = lc['fluxerr'] > 0.
        idx &= lc['flux'] >= 0.
        idx &= lc['snr'] >= 1.

        lc = lc[idx]

        print('after sel')
        print(lc['flux', 'fluxerr', 'zp'])

        if 'sat' in lc.columns:
            print(lc[['band', 'night', 'flux', 'zp',
                      'zpsys', 'fluxerr', 'seeingFwhmEff', 'sat']])

        # trying to fit here
        # outfit = self.fit(lc)

        # print('fit', outfit)

        # printInfo = len(self.SN) > 0
        # self.plot_SN(lcpath, lc, printInfo)

        """
        sigma_z = 1.e-5
        z = lc.meta['z']
        zmeas = z+gauss(0, sigma_z*(1+z))

        self.model.set(z=lc.meta['z'])

        # register bands
        # for visu: only on value per band
        lcb = lc.to_pandas()

        lcb = lcb.groupby(['filter', 'night'])[['airmass', 'pwv',
                                               'ozone', 'aerosol']].median().reset_index()

        lcb['band_cosmo'] = self.telescope.site_name+'::'+lcb['filter']

        self.register_bands_on_the_fly(lcb)

        vv = lc['filter'].tolist()

        vv = list(map(lambda x: self.telescope.site_name+'::' + x, vv))
        lc['band'] = vv
        del lc['filter']

        print('there man', lc['flux', 'fluxerr'])
        result, fitted_model = self.fitIt(lc)
        # print(result)
        # print(fitted_model)
        # fitted_model = None

        """
        lc = lc.to_pandas()
        lc['band'] = lc['band'].str.split('_').str.get(0)

        lc = Table.from_pandas(lc)
        print('kkk', lc['z'])
        
        if fitted_model is not None:
            sncosmo.plot_lc(lc,
                            model=fitted_model,
                            errors=result.errors,
                            zp=25.,
                            yfigsize=9, pulls=False, figtextsize=2.0)
        else:
            sncosmo.plot_lc(lc, xfigsize=9)

        plt.show(block=False)

    def grab_fluxes(self,**pp):
        """
        Function to grab fluxes from SN Ia parameters

        Parameters
        ----------
        pp: dict
            SN parameter dict

        Returns
        -------
        sn_flux : pandas df
            instance of the SNflux class.

        """
        
        from sn_analysis.sn_flux import SNflux
        
        sn_flux = SNflux(pp['x1'][0],pp['color'][0],
                         pp['daymax'][0],pp['z'][0],pp['ebvofMW'][0])
        
        
        return sn_flux
        


    def coadd_lc_deprecated(self, grp,
                 col_means_weighted=[('flux', 'fluxerr')],
                 col_means=['airmass', 'pwv', 'ozone',
                            'aerosol', 'mean_wave', 'zp', 'time', 'snr_m5', 'snr'],
                 col_round=['airmass', 'pwv', 'ozone',
                            'aerosol'],
                 round_vals=[1, 1, 1, 1],
                 col_unique=['zpsys']):
        """
        Method to coadd light-curve points per night/filter

        Parameters
        ----------
        grp : pandas df
            Data to process.
        col_means_weighted : list(str), optional
            list of cols for weighted mean estimation.
            The default is [('flux','fluxerr')].
        col_means : list(str), optional
            list of cols for mean estimation.
            The default is ['airmass','pwv','ozone','aerosol',
                            'mean_wave','zp','time'].
        col_round : list(str), optional
            list of cols to round.
            The default is ['airmass','pwv','ozone',
                            'aerosol','zp','mean_wave'].
        round_vals : list(int), optional
            list of rounding values corresponding to col_round.
            The default is [2,1,1,1,2,2].
        col_unique : list(str), optional
            list of cols with unique value. The default is ['zpsys'].

        Returns
        -------
        astropy table
        output value

        """

        """
        print('in coadd', len(grp))
        print(grp[['flux', 'fluxerr']])
        """

        grp['weight_flux'] = 1./grp['fluxerr']**2

        dictout = {}
        for vv in col_means_weighted:
            pp = vv[0]
            pp_weight = 'weight_{}'.format(pp)
            weight_sum = np.sum(grp[pp_weight])
            mean_weighted = np.sum(grp[pp]*grp[pp_weight])/weight_sum
            dictout[pp] = [mean_weighted]
            dictout[vv[1]] = [1./np.sqrt(weight_sum)]

        for vv in col_means:
            val = grp[vv].mean()
            if vv in col_round:
                idx = col_round.index(vv)
                val = np.round(val, round_vals[idx])

            dictout[vv] = [val]

        for vv in col_unique:
            dictout[vv] = grp[vv].unique().tolist()

        res_df = pd.DataFrame.from_dict(dictout)
        res_df['snr'] = res_df['flux']/res_df['fluxerr']

        """
        print('finally')
        print(res_df[['flux', 'fluxerr']])
        """
        return res_df

    def fitIt_deprecated(self, lc):
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

    def plot_SN_deprecated(self, lcpath, lc, printInfo=False):
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

    def get_text_deprecated(self, meta, ddict=dict(zip(['x1', 'color'], [1, 2]))):
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

    def simple_text_deprecated(self, var, meta, rounding):
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
    def __init__(self, dbDir, dbName, fields, colors, markers, colName):
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
        colors : list(str)
            colors corresponding to fields.
        markers: list(str)
            markers corresponding to fields.
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

        self.ddplot_colors = dict(zip(fields.split(','), colors.split(',')))
        self.ddplot_markers = dict(zip(fields.split(','), markers.split(',')))
        self.colName = colName
        self.dbName = dbName

    def show_stat(self):
        """
        Function to estimate some stats

        Returns
        -------
        res: numpy array
            stat array

        """

        nights = np.unique(self.ddf['night'])
        fields = np.unique(self.ddf['target_name']).tolist()

        rt = []
        for night in nights:
            r = [night]
            idx = self.ddf['night'] == night
            sela = self.ddf[idx]
            for field in fields:
                idxb = sela['target_name'] == field
                selb = sela[idxb]
                r += [len(selb)]
            rt += [r]

        rstat = np.rec.fromrecords(rt, names=['night']+fields)

        # np.sort(rstat, order='night')
        return rstat

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

        if len(sel) == 0:
            print('No observation for this night')
            return

        fig, ax = plt.subplots(figsize=(15, 8))
        figtitle = self.dbName
        figtitle += '\n night {}'.format(night)
        fig.suptitle(figtitle)
        fig.subplots_adjust(right=0.75)

        mjd0 = sel['mjd'].min()
        ttime = (sel['mjd']-mjd0)*24.
        sel = rf.append_fields(sel, 'time', ttime)
        for key, vals in self.ddplot_colors.items():
            idxb = sel[self.colName] == key
            selb = sel[idxb]
            if len(selb) > 0:
                filt = get_filter_alloc(selb)
                label = key.split('DD:')[-1]+('({})'.format(filt))
                ax.plot(selb['time'], selb['filter'], color=vals,
                        linestyle='None', marker=self.ddplot_markers[key],
                        label=label, mfc='None', markersize=10)

        ax.set_xlabel(r'obs time [h]')
        ax.grid(visible=True)
        # ax.legend(bbox_to_anchor=(1., 0.5),
        #          ncol=1, fontsize=12, frameon=False)
        ax.legend(loc='upper left', bbox_to_anchor=(1.0, 0.5),
                  ncol=1, frameon=False, fontsize=15)
        plt.show(block=False)


def get_filter_alloc(data, bands='ugrizy'):
    """
    Function to get the filter allocation

    Parameters
    ----------
    data : numpy array
        Data to process.
    bands : str, optional
        List of filters to consider. The default is 'ugrizy'.

    Returns
    -------
    res : str
        filter allocation.

    """

    r = []
    for b in bands:
        idx = data['filter'] == b
        r.append(str(len(data[idx])))

    res = '/'.join(r)

    return res

class lc_sn:
    def __init__(self,lcDir,lcName,sn_data):
        """
        class to process and plot LC and fluxes from SN params (orig+fitted)

        Parameters
        ----------
        lcDir : str
            LC directory.
        lcName : str
            LC file name.
        sn_data : pandas df
            SN data (orig+fitted params).

        Returns
        -------
        None.

        """
        """
        
        self.ccols = ['x1','color','mb','daymax','z','x0','ebvofMW',
                      'sigma_x1','sigma_color','sigma_mb']
        self.ccols_fit = ['x1_fit','color_fit','mb_fit','t0_fit','z_fit','x0_fit',
                          'ebvofMW','sigma_x1','sigma_color','sigma_mb']
        """
        self.ccols = ['x1','color','daymax','z','x0','ebvofMW',
                     'sigma_x1','sigma_color']
        self.ccols_fit = ['x1_fit','color_fit','t0_fit','z_fit','x0_fit',
                         'ebvofMW','sigma_x1','sigma_color','chisq_red']
        
        self.corresp = dict(zip(self.ccols_fit,self.ccols))
        
        # get lc file

        self.lcs = Read_LightCurve(file_name=lcName, inputDir=lcDir)
        
        # sn_data
        self.sn_data = sn_data
        
    def get_infos(self,lcpath):
        """
        method to get infos

        Parameters
        ----------
        lcpath : str
            lc path to access LC.

        Returns
        -------
        None.

        """
        
        
        lc_plot = {}
        # grab the light curve
        
        lc = self.lcs.get_table(lcpath)
        
        #grab SN fluxes
        pp = {}
        pp_fit={}
        sn_fluxes_fit = Table()
        sn_fluxes_orig = Table()
    
        if len(self.sn_data) > 0:
            self.sn_data['sigma_color'] = np.sqrt(self.sn_data['Cov_colorcolor'])
            #self.sn_data['sigma_mb'] = np.sqrt(self.sn_data['Cov_mbmb'])
            print(self.sn_data.columns.to_list())
            idx = self.sn_data['SNID'] == lcpath
            sn_fit = self.sn_data[idx]
            ppa = sn_fit[self.ccols_fit]
            pp_fit = ppa.to_dict(orient='list')
            ppa = ppa.rename(columns=self.corresp)
            ppf = ppa.to_dict(orient='list')
            sn_flux_fit = self.grab_fluxes(**ppf)
            pp = sn_fit[self.ccols].to_dict(orient='list')
            sn_flux_orig = self.grab_fluxes(**pp)
            
            # get fitted flux here
            sn_fluxes_fit = sn_flux_fit.get_flux(lc)
            sn_fluxes_orig= sn_flux_orig.get_flux(lc)
          
        lc_plot['flux_fit'] = sn_fluxes_fit
        lc_plot['flux_orig'] = sn_fluxes_orig
        
        if pp:
            lc["phase"] = (lc['time']-pp['daymax'][0])/(1.+pp['z'][0])
        lc_plot['lc'] = lc
        
        self.lc_plot = lc_plot
        self.pp_fit = pp_fit
        self.pp = pp
    
    def grab_fluxes(self,**pp):
        """
        method to grab fluxes from SN Ia parameters
    
        Parameters
        ----------
        pp: dict
            SN parameter dict
    
        Returns
        -------
        sn_flux : pandas df
            instance of the SNflux class.
    
        """
        
        from sn_analysis.sn_flux import SNflux
        
        sn_flux = SNflux(pp['x1'][0],pp['color'][0],
                         pp['daymax'][0],pp['z'][0],pp['ebvofMW'][0])
        
        
        return sn_flux
    
    
    def plot_all(self,timescale='phase',fig=None,ax=None):
        """
        method to plot LC+fluxes (original+fit)

        Parameters
        ----------
        timescale : str, optional
            Timescale for the plot. The default is 'phase'.

        Returns
        -------
        None.

        """
        
        #grab the lc to get the filters
        
        lc = self.lc_plot['lc']
        
        """
        idx = lc['fluxerr'] > 0.
        idx &= lc['flux'] >= 0.
        """
        idx = lc['snr'] >= 1.
        
        lc = lc[idx]
        
        bands = np.unique(lc['filter'])
        
        index = dict(zip('ugrizy',[1,2,3,4,5,6]))
        r = []
        for key, vals in index.items():
            r.append((key,vals))
            
        tti = Table(rows=r,names=['filter','index'])
        
        from astropy.table import join
        
        lc = join(lc,tti,keys=['filter'])
        lc['index'] -= np.min(lc['index'])
        
        nfilt_init = len(bands)
        
        if nfilt_init%2==1:
            nfilt = nfilt_init+1
            
        ncols = 2
        nrows = int(nfilt/ncols)
        
        ppos = dict(zip(range(0,6),[(0,0),(0,1),(1,0),(1,1),(2,0),(2,2)]))
        
        if fig is None:
            fig, ax = plt.subplots(nrows=nrows,ncols=ncols,figsize=(11,12))
        
        figtit = ''
        
        if self.pp:
            print(self.pp_fit.keys())
            vinit = ''
            for vv in ['x1','color']:
                x_orig = np.round(self.pp[vv][0],2)
                x_fit = np.round(self.pp_fit['{}_fit'.format(vv)][0],2)
                x_fit_err = np.round(self.pp_fit['sigma_{}'.format(vv)][0],2)
                vinit += '{}={}/{}$\pm$ {}'.format(vv,x_orig,x_fit,x_fit_err)
                if vv == 'x1':
                    vinit += ' - '
            figtit += '{}'.format(vinit)+ os.linesep
            zfit = np.round(self.pp_fit['z_fit'][0],2)
            dfit = np.round(self.pp_fit['t0_fit'][0],2)
            chisq = np.round(self.pp_fit['chisq_red'][0],2)
            t0_str = '$T_0$'
            figtit += 'z={}/{}={}'.format(zfit,t0_str,dfit)+os.linesep
            figtit += '$\chi^2/Ndof$='+'{}'.format(chisq)
            
            
        fig.suptitle(figtit,fontsize=15)
        
        
        index = np.unique(lc['index']).tolist()
        
        sorted(index)
        
        filtercolors = dict(zip('ugrizy', ['b', 'c', 'g', 'y', 'r', 'm']))
        
        for ind in index:
            
            idx = lc['index'] == ind
            """
            idx &= lc['snr'] >= 2
            
            lc['snr_new'] = lc['flux']/lc['fluxerr']
            print('allo',lc[['filter','snr','snr_new']])
            """
            sel_lc = lc[idx]
            b = np.unique(sel_lc['filter'])[0]
            
            ipos = ppos[ind][0]
            jpos = ppos[ind][1]
            
            ax_ = ax[ipos,jpos]
            ax_.errorbar(sel_lc[timescale],
                                   sel_lc['flux'],
                                   yerr=sel_lc['fluxerr'],linestyle='None',
                                   color=filtercolors[b],marker='o',
                                   markersize=5)
            print(sel_lc[['filter','airmass','zp']])
            
            for key, vals in self.lc_plot.items():
                if key != 'lc' and self.pp:
                    idx = vals['filter'] == 'LSST:'+b
                    sel_flux = vals[idx]
                    sel_flux = sel_flux.sort_values(by=[timescale])
                    ls = 'solid'
                    if key == 'flux_orig':
                        ls = 'dotted'
                    ax_.plot(sel_flux[timescale],
                                   sel_flux['flux'],
                                   color=filtercolors[b],linestyle=ls)
            ax_.set_ylabel('flux [pe/s]')      
            ax_.set_xlabel('phase [day]')       
                
            ax_.grid(visible=True)
            
            
            if jpos == 1:
                ax_.yaxis.set_label_position("right")
                ax_.yaxis.tick_right()
            
        #remove empty axes (if any)
        for axr in ax.flat[nfilt_init:]:
            axr.remove()
            
        plt.show(block=False)
        
        
        
        
        
    