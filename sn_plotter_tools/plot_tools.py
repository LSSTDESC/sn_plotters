#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 14:19:57 2026

@author: philippe.gris@clermont.in2p3.fr
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from sn_analysis.sn_tools import get_spline
from . import plt,filtercolors
import pandas as pd
from sn_analysis.sn_tools import fit_lin

def plot_xy(df, 
            xvar='z',xleg='$z$',
            yvar='w',yleg='$w_{DE}$',
            fig=None,ax=None,
            label='',color='k',marker='o',
            linestyle='solid',figtit=''):
    """
    Generic plot function

    Parameters
    ----------
    df : pandas df
        Data to plot.
    xvar : str, optional
        x-axis variable. The default is 'z'.
    xleg : str, optional
        x-axis label. The default is '$z$'.
    yvar : str, optional
        y-axis variable. The default is 'wz'.
    yleg : str, optional
        y-axis label. The default is '$w_{DE}$'.
    fig : matplotlib figure, optional
        figure for the plot. The default is None.
    ax : matplotlib axis, optional
        axis for the plot. The default is None.
    label : str, optional
        plot label. The default is ''.
    color : str, optional
        plot color. The default is 'k'.
    marker : str, optional
        plot marker. The default is 'o'.
    linestyle : str, optional
        plot linestyle. The default is 'solid'.
    figtit : str, optional
        Figure title. The default is ''.

    Returns
    -------
    None.

    """
    
    
    fig_orig=True
    if fig is None:
        fig_orig=False
        fig, ax = plt.subplots(figsize=(12, 8))
        
    if figtit != '':
        fig.suptitle(figtit)
    
    ax.plot(df[xvar], df[yvar], label=label,
               marker=marker, color=color, linestyle=linestyle)
    
    if not fig_orig:
        ax.set_xlabel(r'{}'.format(xleg))
        ax.set_ylabel(r'{}'.format(yleg))
        ax.grid(visible=True)
        ax.legend()
        
def get_grid(tab,varx,vary,varz):
    """
    Function to estimate a grid of data

    Parameters
    ----------
    tab : astropy table
        Data to make the grid from.
    varx : str
        x-axis variable.
    vary : str
        y-axis variable.
    varz : str
        z-axis variable (3rd dimension).

    Returns
    -------
    X : 2D array
        x meshgrid values.
    Y : 2D array
        y meshgrid values.
    Z: : array
        z meshgrid values.

    """
    
    xmin,xmax,xstep,nx = limVals(tab, varx)
    ymin,ymax,ystep,ny = limVals(tab, vary)
    
    xstep = np.round(xstep, 3)
    ystep = np.round(ystep, 3)
    
    xv = np.linspace(xmin, xmax, nx)
    yv = np.linspace(ymin, ymax, ny)
    
    index = np.lexsort((tab[vary], tab[varx]))
    print(tab[index][[varx,vary,varz]])
    flux = np.reshape(tab[index][varz], (nx, ny))
    
    grid=RegularGridInterpolator((xv,yv),flux, 
                                 method='linear', 
                                 bounds_error=False, fill_value=0.)
    #print(grid(()))
    xvp = np.linspace(xmin, xmax, 100*nx)
    yvp = np.linspace(ymin, ymax, 100*ny)
    
    X,Y = np.meshgrid(xvp,yvp)
    
    Z = grid((X,Y))
    
    return X,Y,Z
    
def plot_grid(tab, varx='airmass',xlabel='airmass',
              vary='sigma_pwv',ylabel='$\sigma_{PWV}$ [mm]',
              varz='std_zp_y',figtitle='$\sigma_{ZP}^{y}$',
              iso=[1.,2.,5.],
              txt_iso=['1 mmag','2 mmag','5 mmag'],
              x_iso=[1.5]*3,smoothIt=True):
    """
    Function to make meshgrid plots

    Parameters
    ----------
    tab : astropy table
        Data to process
    varx : str, optional
        x-axis variable. The default is 'airmass'.
    xlabel : str, optional
        x-axis label. The default is 'airmass'.
    vary : str, optional
        y-axis variable. The default is 'sigma_pwv'.
    ylabel : str, optional
        y-axis label. The default is '$\sigma_{PWV}$ [mm]'.
    varz : str, optional
        z-axis variable. The default is 'std_zp_y'.
    figtitle : str, optional
        Figure title. The default is '$\sigma_{ZP}^{y}$'.
    iso : list(float), optional
        List of isocurve variables. The default is [1.,2.,5.].
    txt_iso : str, optional
        List of text for iso curves. The default is ['1 mmag','2 mmag','5 mmag'].
    x_iso : list(float), optional
        x-positions for txt_iso. The default is [1.5]*3.
    smoothIt : bool, optional
        To smooth the iso curves. The default is True.

    Returns
    -------
    None.

    """

    fig,ax = plt.subplots(figsize=(12,8))
    fig.suptitle(figtitle)
    
    #grab the grid
    X,Y,fluxpixels = get_grid(tab,varx,vary,varz)
    
    #grab grid limits
    xmin = np.min(X)
    xmax = np.max(X)
    ymin = np.min(Y)
    ymax = np.max(Y)
    
    #show grid
    im = ax.imshow(fluxpixels,
                   extent=[xmin,xmax,ymin,ymax],
                   #vmin=np.min(fluxpixels),vmax=np.max(fluxpixels),
                   cmap=plt.cm.jet,aspect='auto',origin='lower')
    
    #estimate specific values
    for io,vv in enumerate(iso):
        solutions = np.argwhere((fluxpixels>=vv)&(fluxpixels<=1.1*vv))
        ival = solutions[:,0].tolist()
        jval = solutions[:,1].tolist()
        x_iso = X[ival,jval]
        y_iso = Y[ival,jval]
        df_iso = pd.DataFrame(x_iso,columns=[varx])
        df_iso[vary] = y_iso
        df_iso = df_iso.sort_values(by=[varx])
        df_iso = df_iso.groupby(varx)[vary].mean().reset_index()
      
        
        if not smoothIt:
            ax.plot(df_iso[varx],df_iso[vary],
                    color='k',marker='.',markersize=0.05)
        else:
            if len(df_iso) == 0:
                continue
            xnew, spl_smooth = get_spline(df_iso,varx,vary,nx=10)
            idb = spl_smooth>= ymin
            idb &= spl_smooth <= ymax
            
            
            ax.plot(xnew[idb], spl_smooth[idb],color='k',marker='.',markersize=0.05)
        
            ytext = df_iso[vary].max()+0.00005
            idd = np.argmin(np.abs(df_iso[varx]-x_iso[io]))
            ytext = df_iso.loc[idd,vary]*1.30
            ax.text(1.6,ytext,txt_iso[io])
        
    fig.colorbar(im)
    ax.grid(visible=True)
    ax.set_xlabel(r'{}'.format(xlabel))
    ax.set_ylabel(r'{}'.format(ylabel))
    
def limVals(lc, field):
    """ Get unique values of a field in  a table
    Parameters
    ----------
    lc: Table
        astropy Table (here probably a LC)
    field: str
        name of the field of interest
    
    Returns
    -------
    vmin: float
        min value of the field
    vmax: float
        max value of the field
    vstep: float
        step value for this field (median)
    nvals: int
        number of unique values
    """

    lc.sort(field)
    #dfb = df.sort_values(by=[field])
    vals = np.unique(lc[field].data.round(decimals=4))
    
    vmin = np.min(vals)
    vmax = np.max(vals)
    vstep = np.median(vals[1:]-vals[:-1])

    return vmin, vmax, vstep, len(vals)

def plot_airmass(df,varx='sigma_pwv',xlabel='$\sigma_{PWV}$ [mm]',
                 vary_prefix='std_zp',ylabel='$\sigma_{ZP}$ [mmag]',
                 airmass=[1.2,2.5],
                 y_iso=[1,2,5],
                 txt_iso=['1 mmag','2 mmag','5 mmag'],
                 ymax=6,deltay_txt=0.03,xtext=0.015,
                 smoothIt=False,fitIt=False):
    """
    Function to make a 2D plot for defined airmass values

    Parameters
    ----------
    df : pandas df
        Data to plot.
    varx : str, optional
        x-axis variable. The default is 'sigma_pwv'.
    xlabel : str, optional
        x-axis label. The default is '$\sigma_{PWV}$ [mm]'.
    vary_prefix : str, optional
        y-axis var prefix. The default is 'std_zp'.
    ylabel : str, optional
        y-axis label. The default is '$\sigma_{ZP}$ [mmag]'.
    airmass : list(float), optional
        List of airmass to plot. The default is [1.2,2.5].
    y_iso : list(float), optional
        List of values for isocurves. The default is [1,2,5].
    txt_iso : list(str), optional
        isocurve labels. The default is ['1 mmag','2 mmag','5 mmag'].
    ymax : float, optional
        y-max value. The default is 6.
    deltay_txt : float, optional
        Delta y for the text. The default is 0.03.
    xtext : float, optional
        x value for the text. The default is 0.015.
    smoothIt : bool, optional
        To smooth the data (using spline). The default is False.
    fitIt : bool, optional
        To fit the data (linear). The default is False.    

    Returns
    -------
    None.

    """
    atmos_params = ['sigma_airmass','sigma_aerosol','sigma_pwv','sigma_ozone']
    atmos_paramsb = ['rel_err_airmass','rel_err_aerosol',
                     'rel_err_pwv','rel_err_ozone']
    dx = [0.005,0.001,0.01,5]
    dxb = [0.2,0.5,0.5,1]
    if varx in atmos_params:
        deltax_fit = dict(zip(atmos_params,dx))
    else:
        deltax_fit = dict(zip(atmos_paramsb,dxb))
    
    
    df = df.round({'mean_airmass':2})
    fig, ax = plt.subplots(figsize=(12,8))
    
    bands = 'grizy'
    markers = ['o','P','s','*','h']
    mm = dict(zip(bands,markers))
    lstyle = ['solid','dotted']
    ls = dict(zip(airmass,lstyle))
    
    for airm in airmass:
        idx = df['mean_airmass'] == airm
        sel = df[idx]
        for b in bands:
            yvar = '{}_{}'.format(vary_prefix,b)
            lab = '{} band'.format(b)
            if airm > airmass[0]:
                lab=None
            plot_indiv(sel,xvar=varx,yvar=yvar,label=lab,
                       color=filtercolors[b],
                       marker=mm[b],lstyle=ls[airm],
                       fig=fig,ax=ax,smoothIt=smoothIt,
                       fitIt=fitIt,deltax_fit=deltax_fit[varx])
    
    idx = df['mean_airmass'].isin(airmass)
    sel = df[idx]
    xmin = sel[varx].min()
    xmax = sel[varx].max()
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([0,ymax])
    for io,yvals in enumerate(y_iso):
        ax.plot([xmin,xmax],[yvals]*2,linestyle='dashed',color='k')
        #ax.text(xtext,yvals+deltay_txt,txt_iso[io],fontsize=12)
        #ax.lines(x=xtext, ymin=, ymax=250, color = 'black', linestyles="dashed")
        """
        ax.text(x=xtext, y=yvals, s=txt_iso[io], 
                ha='center', va='center', color='k',
                backgroundcolor='white',fontsize=12)
        """
        ax.text(x=1.05*xmax, y=yvals, s=txt_iso[io], 
                ha='center', va='center', color='k',
                backgroundcolor='white',fontsize=12)
        
    ax.set_xlabel(r'{}'.format(xlabel))
    ax.set_ylabel(r'{}'.format(ylabel))
    ax.legend(loc='upper left',
              bbox_to_anchor=(0., 1.15), ncol=5, frameon=False, fontsize=15)
    ax.grid(visible=True)
    #ax.text(0.2,1.05,'.... airmass=2.5 ',fontsize=12,transform=ax.transAxes)
    
    x_trans=0.25
    ax.annotate('', xy=(x_trans+0.,1.05), 
                xycoords='axes fraction', xytext=(x_trans+0.05, 1.05),
                arrowprops=dict(arrowstyle="-", color='k'))
    ax.text(x_trans+0.055,1.04,'airmass={}'.format(airmass[0]),
            fontsize=12,transform=ax.transAxes)
    ax.annotate('', xy=(x_trans+0.2,1.05), xycoords='axes fraction',
                xytext=(x_trans+0.25, 1.05),
               arrowprops=dict(arrowstyle="-", color='k',linestyle='dotted'))
    ax.text(x_trans+0.255,1.04,'airmass={}'.format(airmass[1]),
            fontsize=12,transform=ax.transAxes)
    
def plot_indiv(df,
               xvar='z', xlabel='z', 
               yvar='N', ylabel='NSN',label='',
               lstyle='solid',color='k',marker='o',
               figtitle='',fig=None, ax=None,
               smoothIt=False,fitIt=False,deltax_fit=0.01):
    """
    Function to make indiv plot

    Parameters
    ----------
    df : pandas df
        Data to plot.
    xvar : str, optional
        x-axis variable. The default is 'z'.
    xlabel : str, optional
        x-axis label. The default is 'z'.
    yvar : str, optional
        y-axis variable. The default is 'N'.
    ylabel : str, optional
        y-axis label. The default is 'NSN'.
    label : str, optional
        label. The default is ''.
    lstyle : str, optional
        Line style. The default is 'solid'.
    color : str, optional
        color. The default is 'k'.
    marker : str, optional
        marker. The default is 'o'.
    figtitle : str, optional
        figure title. The default is ''.
    fig : matplotlib figure, optional
        Figure for the plot. The default is None.
    ax : matplotlib axis, optional
        Axis for the plot. The default is None.
    smoothIt : bool, optional
        To smooth the plot using spline. The default is False.
    fitIt : bool, optional
        To fit the plot using linear function. The default is False.
    Returns
    -------
    None.

    """
    
    if fig is None:
        fig, ax = plt.subplots(figsize=(12,8))
    
    
    if not smoothIt and not fitIt:
        ax.plot(df[xvar],df[yvar],
                marker=marker,linestyle=lstyle,
                color=color,markersize=8,mfc='None',label=label)
    if smoothIt:
        df = df.sort_values(by=[xvar])
        xnew, spl_smooth = get_spline(df,xvar,yvar,nx=10)
        ax.plot(xnew, spl_smooth,color=color,
                marker=marker,linestyle=lstyle,
                markersize=10,mfc='None',label=label) 
        
    if fitIt:
        res = fit_lin(df,xvar,yvar)
        xmin=df[xvar].min()
        xmax=df[xvar].max()
        xv = np.arange(xmin,xmax+deltax_fit,deltax_fit)
        yv = res[0]*xv+res[1]
        ax.plot(xv,yv,color=color,
                marker=marker,linestyle=lstyle,
                markersize=10,mfc='None',label=label,lw=2)