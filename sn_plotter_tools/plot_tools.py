#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 14:19:57 2026

@author: philippe.gris@clermont.in2p3.fr
"""
import matplotlib.pyplot as plt

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
