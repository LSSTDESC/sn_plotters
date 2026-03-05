#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 13:19:36 2026

@author: philippe.gris@clermont.in2p3.fr
"""
import sncosmo
from sn_telmodel.sn_throughputs import get_telescope
import numpy as np
from astropy.table import Table,vstack

class SNflux:
    def __init__(self,x1,color,x0,daymax,z,
                 model='salt3',
                 version='2.0',
                 absmag=-19.0906,
                 magsys='vega',band='bessellB',
                 tel_dir = 'throughputs',
                 throughputsDir = 'baseline',
                 atmosDir = 'atmos',
                 tag_tel = '1.9',
                 airmass=1.2,
                 aerosol=0.0,
                 pwv=4.0,
                 ozone=400,):
        """
        class to estimate SN Ia  flux vs time

        Parameters
        ----------
        x1 : float
            SN stretch parameter.
        color : float
            SN color parameter.
        x0 : float
            SN x0 parameter.
        daymax :  float
            SN daymax (T0) parameter.
        z : float
            SN redshift parameter.
        model : str, optional
            sncosmo model. The default is 'salt3'.
        version : float, optional
            sncosmo model version. The default is 2.0.
        absmag : float, optional
            SN abs mag. The default is -19.0906.
        magsys : str, optional
            mag system. The default is 'vega'.
        band : str, optional
            SN band. The default is 'bessellB'.
        tel_dir : str, optional
            telescope directory. The default is 'throughputs'.
        throughputsDir : str, optional
            throughput directory. The default is 'baseline'.
        atmosDir : str, optional
            atmosphere directory. The default is 'atmos'.
        tag_tel : str, optional
            telescope tag version. The default is '1.9'.
        airmass : float, optional
            airmass default value. The default is 1.2.
        aerosol : float, optional
            aerosol default value. The default is 0.0.
        pwv : float, optional
            pwv default value. The default is 4.0.
        ozone : float, optional
            ozone default value. The default is 400.

        Returns
        -------
        None.

        """
        
        
        # SN parameters
        self.x1 = x1
        self.color = color
        self.x0 = x0
        self.daymax = daymax
        self.z = z
        self.model = model
        self.version = version
        self.absmag = absmag
        self.band = band
        self.magsys = magsys
        
        
        #telescope parameters
        self.tel_dir = tel_dir
        self.throughputsDir = throughputsDir
        self.atmosDir = atmosDir
        self.tag_tel = tag_tel
        
        #atmos parameters
        self.airmass = airmass
        self.pwv = pwv
        self.aerosol = aerosol
        self.ozone = ozone
        
        #instances of SN and telescope
        
        self.sn = self.get_sn()
        self.telescope = self.get_telescope()
        
    def get_sn(self):
        """
        Method for SN Ia instance in sncosmo

        Returns
        -------
        sn : sncosmo.Model
            SN Ia instance.

        """
        
        source = sncosmo.get_source(self.model, self.version)
        dustmap = sncosmo.OD94Dust()
        sn = sncosmo.Model(source=source,
                           effects=[dustmap, dustmap],
                           effect_names=['host', 'mw'],
                           effect_frames=['rest', 'obs'])
        sn.set(z=self.z)
        sn.set(t0=self.daymax)
        sn.set(x1=self.x1)
        sn.set(c=self.color)
        sn.set(x0=self.x0)
        
        return sn
    
    def get_telescope(self):
        """
        Method for a telescope instance

        Returns
        -------
        telescope : Throughput class
            Total Throughput.

        """
        
        telb = '{}_{}'.format(self.tel_dir, self.tag_tel)
        through_dir = '{}/{}'.format(telb, self.throughputsDir)
        atmos_dir = '{}/{}'.format(telb, self.atmosDir)
        telescope = get_telescope(tel_dir=telb,
                                  through_dir=through_dir,
                                  atmos_dir=atmos_dir,
                                  tag=self.tag_tel, 
                                  airmass=self.airmass,
                                  aerosol=self.aerosol, 
                                  pwv=self.pwv, 
                                  ozone=self.ozone)
         
        return telescope
        
        
    def get_flux(self,lc_data):
        """
        Method to get sn_fluxes from LC (time+filters)

        Parameters
        ----------
        lc_data : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        
        print('oooooo',lc_data.columns)
        tmin = self.daymax-20*(1+self.z)
        tmax = self.daymax+60*(1+self.z)
        
        ccols = ['time','band_cosmo','filter','airmass','pwv','aerosol','ozone']
        lc = lc_data[ccols]
        
        filters = np.unique(lc['filter'])
        
        airmassb = np.round(self.airmass,2)
        pwvb = np.round(self.pwv,2)
        aerosolb = np.round(self.aerosol,2)
        ozoneb = np.round(self.ozone,2)
        
        b_filt = {}
        for fi in filters:
            bcols = self.telescope.site_name+'::'
            bcols += fi+'_'
            bcols += '{}'.format(airmassb)+'_' 
            bcols +='{}'.format(pwvb)+'_'      
            bcols += '{}'.format(ozoneb)+'_' 
            bcols +='{}'.format(aerosolb)
            b_filt[fi] = bcols
            
        r = []
         
        tis = np.arange(tmin,tmax,0.5)
        
        import pandas as pd
        df_add = pd.DataFrame()
        for key,vals in b_filt.items():
            
            df = pd.DataFrame(tis, columns=['time'])
            df['band_cosmo'] = vals
            df['filter'] = fi
            df['airmass'] = airmassb
            df['pwv'] = pwvb
            df['ozone'] = ozoneb
            df['aerosol'] = aerosolb
            df_add = pd.concat((df_add,df))
          
        print(df_add)
        
        lc_tot = vstack([lc,Table.from_pandas(df_add)])
        
        print(lc_tot)
    def register_bands(self,data):
        from sn_tools.sn_utils import register_bands_sncosmo
        
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
        
        
        
        
        
        