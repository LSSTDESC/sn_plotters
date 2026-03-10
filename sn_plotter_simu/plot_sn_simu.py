#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 08:47:26 2026

@author: philippe.gris@clermont.in2p3.fr
"""
import matplotlib.pyplot as plt

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
    print(sn_flux)
    mjd_min = sn_flux['phase'].min()
    mjd_max = sn_flux['phase'].max()
    flux_min={}
    flux_max={}
    
    for b in bands:
        idx = sn_flux['filter'] == 'LSST:'+b
        sel = sn_flux[idx]
        flux_min[b] = sel['flux'].min()
        flux_max[b] = sel['flux'].max()
        
    for sed in sn_sed:
        fig = plt.figure(figsize=(12,8))
        ax1 = fig.add_subplot(2,1,1)
        mjd = sed.meta['mjd']
        fig.suptitle('MJD:{}'.format(mjd))
        ax1.plot(sed['wavelength'],sed['flux'],'k.')
        ax1.grid(visible=True)
        ax1.set_xlim([4500.,20000.])
        for i,b in enumerate(bands):
            idx = sn_flux['filter'] == 'LSST:'+b
            idx &= sn_flux['time'] <= mjd
            sel_flux = sn_flux[idx]
            ax = fig.add_subplot(2,3,i+4)
            ax.plot(sel_flux['phase'],sel_flux['flux'],'ko')
            ax.set_xlim([mjd_min,mjd_max])
            ax.set_ylim([flux_min[b],flux_max[b]])
            ax.grid(visible=True)
        plt.show()

