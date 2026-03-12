#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 08:47:26 2026

@author: philippe.gris@clermont.in2p3.fr
"""
#import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from . import plt,filtercolors

def plot_flux_spectra(sn_flux,sn_sed,outDir='None',bands='izy',phase_to_draw=[]):
    """
    Function to plot spectra (top) and flux from SN Ia.

    Parameters
    ----------
    sn_flux_o : astropy table
        SN flux.
    sn_sed : list of astropytables
        SN spectra.
    outDir: str, optional.
        Dir where to save data. The default is None

    Returns
    -------
    None.

    """
    if not phase_to_draw:
        phases = np.unique(sn_sed['phase']).tolist()
    else:
        phases = get_phase(sn_sed,phase_to_draw)
    
    #sed_max = sn_sed['flux'].max()
    
    if phase_to_draw:
        fig, ax = plt.subplots(nrows=2,figsize=(12,8))
    
    figtitm = '(x1,color,z)=({},{},{}) \n'.format(sn_flux.meta['x1'],
                                              sn_flux.meta['color'],
                                              sn_flux.meta['z'])
    for ip,phase in enumerate(phases):
        
        figtit = figtitm
        #fig = plt.figure(figsize=(12,8))
        if not phase_to_draw:
            fig, ax = plt.subplots(nrows=2,figsize=(12,8))
        #ax1 = fig.add_subplot(2,1,1)
        mjd = plot_sed(sn_sed,phase,figtit,fig=fig,ax=ax[0])
        #ax1.set_xlim([0.,sed_max])
        mjda=9.*10**9
        label = True
        if not phase_to_draw:
            mjda = mjd
        else:
            if ip > 0:
                continue
        
        plot_flux(sn_flux,bands,mjda,fig=fig,ax=ax[1],labelIt=label)
        
        
        if outDir != 'None': 
            fName = '{}/flux_spectra_{}.png'.format(outDir,str(ip).zfill(3))
            plt.savefig(fName)
            plt.close(fig)

    if phase_to_draw:
        plt.show()

def plot_sed(sn_sed,phase,figtit,fig=None,ax=None):
    
    idx = sn_sed['phase'] == phase
    sel_sed = Table(sn_sed[idx])
    mjd = np.unique(sel_sed['mjd'])[0]
    mjd_str = 'MJD:{}'.format(mjd)
    figtit += '{}'.format(mjd_str)
    fig.suptitle(figtit,fontweight='bold')
    ax.plot(sel_sed['wavelength']/10.,sel_sed['flux'],
             color='k',marker='.',markersize=3,label='phase={}'.format(phase))
    ax.grid(visible=True)
    ax.set_xlim([450.,1300.])
    ax.set_xlabel('wavelength [nm]')
    ax.set_ylabel('SED [erg/s/cm$^2$/A]')
    ax.legend()
    
    return mjd
    
def plot_flux(sn_flux,bands,mjd=9.10**9,
              varx='phase',legx='phase',
              fig=None,ax=None,labelIt=True):
    
    
    if fig is None:
        fig, ax = plt.subplots(figsize=(12,8))
        
    idx = np.in1d(sn_flux['filter_notel'],list(bands))
    sn_flux = sn_flux[idx]
    
    mjd_min = sn_flux[varx].min()
    mjd_max = sn_flux[varx].max()
    
    flux_min = sn_flux['flux'].min()
    flux_max = sn_flux['flux'].max()
    
    
    for i,b in enumerate(bands):
        idx = sn_flux['filter'] == 'LSST:'+b
        idx &= sn_flux['time'] <= mjd
        sel_flux = sn_flux[idx]
        #ax = fig.add_subplot(2,3,i+4)
        #ax = fig.add_subplot(2,1,2)
        thelabel = b
        if not labelIt:
            thelabel=None
        ax.plot(sel_flux['phase'],sel_flux['flux'],
                 color=filtercolors[b],label=thelabel)
        ax.set_xlim([mjd_min,mjd_max])
        #ax2.set_ylim([flux_min[b],flux_max[b]])
        ax.set_ylim([flux_min,1.05*flux_max])
        ax.grid(visible=True)
        ax.legend()
    ax.set_xlabel(legx)
    ax.set_ylabel('flux [pe/s]')
 
def get_phase(df,phase_list):
    
    r = []
    for ph in phase_list:
        df['diff'] = np.abs(df['phase']-ph)
        idx = np.argmin(df['diff'])
        print('ohohoh',ph,df[idx]['phase'])
        r.append(df[idx]['phase'])
        
    return r
    
    