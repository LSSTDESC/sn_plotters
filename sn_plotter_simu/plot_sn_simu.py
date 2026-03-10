#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 08:47:26 2026

@author: philippe.gris@clermont.in2p3.fr
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

def plot_flux_spectra(sn_flux,sn_sed):
    """
    Function to plot spectra (top) and flux from SN Ia.

    Parameters
    ----------
    sn_flux : astropy table
        SN flux.
    sn_sed : list of astropytables
        SN spectra.

    Returns
    -------
    None.

    """
    
    bands = 'izy'
    
    """
    idx = sn_flux['filter'].isin(bands)
    sel_flux = sn_flux[idx]
    """
    
    mjd_min = sn_flux['phase'].min()
    mjd_max = sn_flux['phase'].max()
    """
    flux_min={}
    flux_max={}
    
    for b in bands:
        idx = sn_flux['filter'] == 'LSST:'+b
        sel = sn_flux[idx]
        flux_min[b] = sel['flux'].min()
        flux_max[b] = sel['flux'].max()
    """
    flux_min = sn_flux['flux'].min()
    flux_max = sn_flux['flux'].max()
    print('sed plotting',len(sn_sed))
    phases = np.unique(sn_sed['phase']).tolist()
    print(phases)
    
    for phase in phases:
        idx = sn_sed['phase'] == phase
        sel_sed = Table(sn_sed[idx])
        #fig = plt.figure(figsize=(12,8))
        fig, ax = plt.subplots(nrows=2,figsize=(12,8))
        #ax1 = fig.add_subplot(2,1,1)
        ax1= ax[0]
        mjd = np.unique(sel_sed['mjd'])[0]
        fig.suptitle('MJD:{}'.format(mjd))
        ax1.plot(sel_sed['wavelength'],sel_sed['flux'],'k.')
        ax1.grid(visible=True)
        ax1.set_xlim([4500.,20000.])
        ax1.set_xlabel('wavelength []')
        ax1.set_ylabel('flux []')
        ax2 = ax[1]
        for i,b in enumerate(bands):
            idx = sn_flux['filter'] == 'LSST:'+b
            idx &= sn_flux['time'] <= mjd
            sel_flux = sn_flux[idx]
            #ax = fig.add_subplot(2,3,i+4)
            #ax = fig.add_subplot(2,1,2)
            ax2.plot(sel_flux['phase'],sel_flux['flux'],'ko')
            ax2.set_xlim([mjd_min,mjd_max])
            #ax2.set_ylim([flux_min[b],flux_max[b]])
            ax2.set_ylim([flux_min,1.1*flux_max])
            ax2.grid(visible=True)
        ax2.set_xlabel('phase')
        ax2.set_ylabel('flux [pe/s]')
        plt.show()

