#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 16:07:43 2026

@author: philippe.gris@clermont.in2p3.fr
"""
from . import plt
from sn_plotter_tools.plot_tools import plot_airmass
import pandas as pd

atmos_params = ['airmass','pwv','ozone','aerosol','total']
markers = ['o','s','P','h','v']
marks = dict(zip(atmos_params,markers))
    
def plot_summary_obs_param(df,obs_param='zp',airmass=1.2,
         fig=None,ax=None,labelIt=True,
         linestyle='solid',color='k',valref=1):
    """
    Function to plot zp,mean_wave vs inverse of the fit slope

    Parameters
    ----------
    df : pandas df
        Data to plot.
    obs_param : str, optional
        obs parameter. The default is 'zp'.
    airmass : float, optional
        airmass value. The default is 1.2.
    fig : matplotlib figure, optional
        Figure for the plot. The default is None.
    ax : matplotlib axis, optional
        axis for the plot. The default is None.
    labelIt : bool, optional
        To add a label to the plot. The default is True.
    linestyle : str, optional
        line style for the plot. The default is 'solid'.
    color : str, optional
        color for the plot. The default is 'k'.
    valref: float, optional.
        reference value. The default is 1.

    Returns
    -------
    None.

    """
    
    
    if fig is None:
        fig,ax = plt.subplots(figsize=(12,8))
    
    idx = df['obs_param'] == obs_param
    idx &= df['airmass'] == airmass
    
    sel = df[idx]
    
    atmos_params = sel['atmos_param'].unique()
    
    #fig, ax = plt.subplots(figsize=(12,8))
    #fig.suptitle('airmass={}'.format(airmass))
    
    atm_ref = ['airmass','pwv','ozone','aerosol']
    markers = ['o','s','P','h']
    marks = dict(zip(atm_ref,markers))
    
    for atm in atmos_params:
        #fig, ax = plt.subplots(figsize=(12,8))
        label = None
        if labelIt:
            label = atm
        idx = sel['atmos_param'] == atm
        selb = sel[idx]
        ax.plot(selb['band'],valref/selb['slope'],
                label=label,linestyle=linestyle,
                marker=marks[atm],mfc='None',markersize=12,color=color)
        
def plot_summary(df_zp,obs_param='zp',valref=1,unit='mmag'):
    """
    summary plot

    Parameters
    ----------
    df_zp : pandas df
        Data to plot.
    obs_param : str, optional
        obs parameter. The default is 'zp'.
    valref : float, optional
        reference value. The default is 1.
    unit : str, optional
        plot unit. The default is 'mmag'.

    Returns
    -------
    None.

    """
    
    fig, ax = plt.subplots(figsize=(12,8))  
    plot_summary_obs_param(df_zp,obs_param=obs_param,fig=fig,ax=ax,valref=valref)
    plot_summary_obs_param(df_zp,obs_param=obs_param,airmass=2.0,
         fig=fig,ax=ax,labelIt=False,linestyle='dotted',color='b',
         valref=valref)
    
    ax.legend(loc='upper left',
                  bbox_to_anchor=(0.1, 1.15), ncol=5, frameon=False, fontsize=15)
    
    ax.grid(visible=True)
    
    ax.set_yscale("log")
    x_trans=0.25
    ax.annotate('', xy=(x_trans+0.,1.05), 
                xycoords='axes fraction', xytext=(x_trans+0.05, 1.05),
                arrowprops=dict(arrowstyle="-", color='k'))
    ax.text(x_trans+0.055,1.04,'airmass=1.2',
            fontsize=12,transform=ax.transAxes)
    ax.annotate('', xy=(x_trans+0.2,1.05), xycoords='axes fraction',
                xytext=(x_trans+0.25, 1.05),
               arrowprops=dict(arrowstyle="-", color='k',linestyle='dotted'))
    ax.text(x_trans+0.255,1.04,'airmass=2.0',
            fontsize=12,transform=ax.transAxes)
    
    ax.set_xlabel(r'band')
    ylabel = '$\sigma_{atmos\ param}$'
    po = obs_param.replace('_','\ ')
    
    ylabel += '($\sigma_{'+po+'}$='+'{}'.format(valref)+' {}'.format(unit)+')'
    ax.set_ylabel(r'{}'.format(ylabel))
    
    #add auxtel performance
    xmin, xmax = ax.get_xlim()
    yv_auxtel=[0.2,20,3e-3,5e-3]
    coeff= [1.10]*2+[0.7]+[1.15]
    ttxt = ['$\sigma_{PWV}$','$\sigma_{ozone}$',
            '$\sigma_{airmass}$','$\sigma_{aerosol}$']
    units = ['mm','DU','','']
    for io,yy in enumerate(yv_auxtel):
        r = []
        for b in 'grizy':
            r.append((b,yy))
            
        dfaux = pd.DataFrame(r,columns=['band','auxres'])
        
        ax.plot(dfaux['band'],dfaux['auxres'],color='r',linestyle='dashed')
        """
        ax.text(x_trans+0.055,ypos[io],ttxt[io]+'='+'{}'.format(yy),
            fontsize=12,transform=ax.transAxes)
        """
        ax.text(3.2,coeff[io]*yy,ttxt[io]+'='+'{}'.format(yy)+' {}'.format(units[io]),
            fontsize=12,color='r')    

def plot_all_summary(df_zp,df_wave):
    """
    plot all summary

    Parameters
    ----------
    df_zp : pandas df
        zp data.
    df_wave : pandas df
        mean wave data.

    Returns
    -------
    None.

    """
    
    plot_summary(df_zp)
    plot_summary(df_wave,obs_param='mean_wave',valref=0.1,unit='nm')
    
def plot_atmos_data_airmass(theDir,
                            atmos_params=['pwv','aerosol','airmass','ozone']):
    """
    plot atmos data

    Parameters
    ----------
    theDir : str
        Data dir.
    atmos_params : list(str), optional
        List of atmos parameters. 
        The default is ['pwv','aerosol','airmass','ozone'].

    Returns
    -------
    None.

    """
    
    all_atm = ['pwv','aerosol','airmass','ozone']
    legxx = ['$\sigma_{PWV}$ [mm]',
            '$\sigma_{aerosol}$',
            '$\sigma_{airmass}$',
            '$\sigma_{ozone}$ [DU]']
    legxxrel = ['$\\frac{\sigma_{PWV}}{<PWV>}$ [%]',
            '$\\frac{\sigma_{aerosol}}{<aerosol>}$ [%]',
            '$\\frac{\sigma_{airmass}}{<airmass>}$ [%]',
            '$\\frac{\sigma_{ozone}}{<ozone>}$ [%]']
    
    airmass=[1.2,2.0]
    xt = [0.15,0.011,0.01,25]
    xxtext=dict(zip(all_atm,xt))
    legx = dict(zip(all_atm,legxx))
    legxrel = dict(zip(all_atm,legxxrel))
    
    for vv in atmos_params:
        theFile = 'zp_atmos_{}.hdf5'.format(vv)
        fName = '{}/{}'.format(theDir,theFile)

        df = pd.read_hdf(fName)

        if vv == 'aerosol':
            idx = df['sigma_aerosol'] <= 0.0125
            df = df[idx]
            
        if vv == 'airmass':
           idx = df['sigma_airmass'] <= 0.04
           df = df[idx]   


        for b in 'grizy':
            df['std_zp_{}'.format(b)] *= 1000 # in mmag
         
        
        
        plot_airmass(df,varx='sigma_{}'.format(vv),xlabel=legx[vv],
                         vary_prefix='std_zp',airmass=airmass, 
                         y_iso=[1,2,3,5],
                         txt_iso=['1 mmag','2 mmag','3 mmag','5 mmag'],
                         xtext=xxtext[vv],smoothIt=False,fitIt=True) 
        
        plot_airmass(df,varx='sigma_{}'.format(vv),xlabel=legx[vv],
                         vary_prefix='std_mean_wave',
                         ylabel='$\sigma_{meanwave}$ [mm]',
                         airmass=airmass,
                         y_iso=[0.05,0.1,0.15],
                         txt_iso=['0.05 nm','0.1 nm','0.15 nm'],
                         ymax=0.2,deltay_txt=0.005,
                         xtext=xxtext[vv],smoothIt=False,fitIt=True)
        
        rel_err = 'rel_err_{}'.format(vv)
        df[rel_err] = 100.*df['sigma_{}'.format(vv)]/df['mean_{}'.format(vv)]
        
        print(df.columns,legxrel[vv])
        plot_airmass(df,varx=rel_err,xlabel=legxrel[vv],
                         vary_prefix='std_zp',airmass=airmass, 
                         y_iso=[1,2,3,5],
                         txt_iso=['1 mmag','2 mmag','3 mmag','5 mmag'],
                         xtext=xxtext[vv],smoothIt=False,fitIt=True)
        plot_airmass(df,varx=rel_err,xlabel=legxrel[vv],
                         vary_prefix='std_mean_wave',
                         ylabel='$\sigma_{meanwave}$ [mm]',
                         airmass=airmass,
                         y_iso=[0.05,0.1,0.15],
                         txt_iso=['0.05 nm','0.1 nm','0.15 nm'],
                         ymax=0.2,deltay_txt=0.005,
                         xtext=xxtext[vv],smoothIt=False,fitIt=True)
        
        
def plot_perf(data,x_main='sigma',
              obs_param='zp',unit='mmag',ylabel='zp',
              atmos_params=['airmass','ozone','aerosol','pwv','total'],
              airmass=[1.2,2.],ylines=[1,5,10],
              yannot=['1 mmag','5 mmag','10 mmag'],extra_leg=''):
    """
    Function to draw a perf plot

    Parameters
    ----------
    data : pandas df
        Data to plot.
    x_main : str, optional
        prefix var. The default is 'sigma'.
    obs_param : str, optional
        obs param. The default is 'zp'.
    unit : str, optional
        obs param unit. The default is 'mmag'.
    ylabel : str, optional
        y-axis label. The default is 'zp'.
    atmos_params : list(str), optional
        list of atmos params. 
        The default is ['airmass','ozone','aerosol','pwv','total'].
    airmass : list(float), optional
        List of airmass to consider. The default is [1.2,2.].
    ylines : list(float), optional
        y-values of lines to draw. The default is [1,5,10].
    yannot : list(str), optional
        ylines annot. The default is ['1 mmag','5 mmag','10 mmag'].
    extra_leg: str, optional.
        extra legend to add to the plot. The default is ''.

    Returns
    -------
    None.

    """
    
    fig, ax = plt.subplots(figsize=(12,8))
    
    ls = dict(zip(airmass,['solid','dotted']))
    color = dict(zip(airmass,['black','red']))
    
    for airm in airmass:
        idx = data['airmass'] == airm
        sel = data[idx]
        
        for atm_param in atmos_params:
            label =atm_param
            if airm > 1.2:
                label = None
            obs_str = '{}_{}_{}'.format(x_main,obs_param,atm_param)
            ax.plot(sel['band'],sel[obs_str],
                    marker=marks[atm_param],mfc='None',
                    linestyle=ls[airm],color=color[airm],label=label)
    
    xmin,xmax = ax.get_xlim()
    for io,yl in enumerate(ylines):
        ax.plot([xmin,xmax],[yl]*2,linestyle='dashed',color='b')
        ax.text(1.02*xmax,yl,'{}'.format(yannot[io]),
            fontsize=12,color='b')       
    ax.grid(visible=True)
    
    if x_main == 'sigma':
        ylabel = '$\\'+x_main+'_{'+ylabel+'}$ ['+unit+']'
    else:
        #ylabel = x_main+'$_{\sigma_{'+ylabel+'}}$ ['+unit+']'
        ylabel = '$\sigma_{'+ylabel+'}$'+' budget ['+unit+']'
    ax.set_ylabel(r'{}'.format(ylabel))
    ax.set_xlabel(r'band')
    
    ax.legend(loc='upper left',
                 bbox_to_anchor=(0.1, 1.15), 
                 ncol=5, frameon=False, fontsize=15)
    
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
    
    ax.set_xlim([xmin,xmax])
    
    if extra_leg != '':
       ax.text(-0.15,0.97,extra_leg,
            fontsize=12,transform=ax.transAxes,color='b') 
    
def plot_perf_obs_param(df_zp,sigmas,unit_atmos,
                    obs_param='zp',unit='nm',ylabel='zp',
                    ylines=[1,5,10],
                    yannot=['1 mmag','5 mmag','10 mmag']):
    """
    Function to draw perf plots for an obs_param

    Parameters
    ----------
    df_zp : pandas df
        Data to plot.
    sigmas : dict
        sigmas of atmos params.
    unit_atmos : dict
        unit for sigmas of atmos params.    
    obs_param : str, optional
        obs param. The default is 'zp'.
    unit : str, optional
        obs param unit. The default is 'nm'.
    ylabel : str, optional
        y-axis label. The default is 'zp'.
    ylines : list(float), optional
        y-values of lines to draw. The default is [1,5,10].
    yannot : list(str), optional
        y-lines annot. The default is ['1 mmag','5 mmag','10 mmag'].

    Returns
    -------
    None.

    """
    
    bands = 'grizy'
    b_index = [0,1,2,3,4]
    dfb = pd.DataFrame(list(bands),columns=['band'])
    dfb['band_index'] = b_index
  
    
    extra_leg = ''
    
    for key,vals in sigmas.items():
        vvar = '$\sigma_{'+key+'}$='+'{}'.format(vals)
        extra_leg += vvar + ' ' +unit_atmos[key]+'\n'
        
    from sn_analysis.sn_atmos_tools import get_atmos_data,get_values
    rr_zp = get_atmos_data(df_zp)
    res_zp = get_values(rr_zp,sigmas=sigmas)
  
    res_zp = res_zp.merge(dfb,left_on=['band'],right_on=['band'])
  
    res_zp = res_zp.sort_values(by=['band_index'])
  
  
    res_zp['frac_check'] = 0
    vart = 'sigma_{}_total'.format(obs_param)
    for atm_param in atmos_params:
        fracx = 'frac_{}_{}'.format(obs_param,atm_param)
        varx = 'sigma_{}_{}'.format(obs_param,atm_param)
        res_zp[fracx] = 100.*res_zp[varx]**2/res_zp[vart]**2
        res_zp['frac_check'] += res_zp[fracx]
  
    plot_perf(res_zp,obs_param=obs_param,unit =unit,
            ylabel=ylabel,ylines=ylines,yannot=yannot,extra_leg=extra_leg)
    
    plot_perf(res_zp,x_main='frac',obs_param=obs_param,unit='%',ylabel=ylabel,
              atmos_params=['airmass','ozone','aerosol','pwv'],ylines=[],
              extra_leg=extra_leg)
    
def plot_results_config(df,config_df,obs_param='zp',unit='mmag',plotDir='',
                 tagline=[1,2,5,10]):
    """
    Function to plot the results

    Parameters
    ----------
    df : pandas df
        Data to plot.
    obs_param : str, optional
        obs parameter to plot (zp/mean_wave). The default is 'zp'.
    unit : str, optional
        unit corresponding to obs_param (mmag/nm). The default is 'mmag'.
    plotDir : str, optional
        Output dir for the plots. The default is ''.
    tagline : list(float), optional
        Lines to add to the plot. The default is [1,2,5,10].

    Returns
    -------
    None.

    """

    from sn_analysis.sn_atmos_tools import add_index_band,add_legend,get_str
    #add index
    df = add_index_band(df)
    
    fig, ax = plt.subplots(figsize=(12,8))
    fig.subplots_adjust(top=0.85,right=0.9)
    configs = df['config'].unique()
    airmass = df['airmass'].unique().tolist()
    airmass = list(map(float, airmass))
    
    yvar= 'sigma_{}_tot'.format(obs_param)
    
    lstyles = dict(zip(airmass,['solid','dotted']))
    markers = ['o','s','P','h','v']
    
    mmarks = dict(zip(configs,markers[:len(configs)]))
    
    for airm in airmass:
        idx = df['airmass'] == airm
        sel = df[idx]
        for config in configs:
            idxb = sel['config'] == config
            selb = sel[idxb]
            idxb = config_df['config']==config
            sel_config = config_df[idxb]
            label = get_str(sel_config,atmos_params=['airmass','ozone',
                                                     'aerosol','pwv'])
            label = '$\sigma_{atmos}$='+label
            if airm > 1.5:
                label = None
            ax.plot(selb['band'],selb[yvar],
                    marker=mmarks[config],mfc='None',
                    linestyle=lstyles[airm],label=label)
            
    ax.grid(visible=True)
    
    ax.legend(loc='upper left',
              bbox_to_anchor=(0.05, 1.2), ncol=2, frameon=False, fontsize=15)
    ax.set_xlabel(r'band')
    ylabel = '$\sigma_{'+obs_param+'}$'+ '[{}]'.format(unit)
    ax.set_ylabel(r'{}'.format(ylabel))
    xmin,xmax = ax.get_xlim()
    for tt in tagline:
        ax.plot([xmin,xmax],[tt]*2,linestyle='dashed',color='k')
        ax.text(1.02*xmax,tt,'{}'.format(tt)+' mmag',fontsize=12)
    ax.set_xlim([xmin,xmax])
    
    # add airmass legend
    add_legend(ax, airmass)
    
    if plotDir != '':
        outName = '{}/summary_zp.png'.format(plotDir)
        plt.savefig(outName)
    
def plot_sigma_obs_param(res, obs_param='zp',unit_obs_param='mmag',
                         vary='sigma_atmos_param',err_rel=False,
                         atmos_params = ['airmass','ozone','aerosol','pwv'],
                         limy =[[0.,0.05],[0.,100.],[0.,0.0055],[0.,0.5]],
                         auxtel_data=[3.e-3,20,5.e-3,0.2],plotDir=''):
    """
    Plot sigma_atmos vs band for a set of sigma_zp_atmos values

    Parameters
    ----------
    res : pandas df
        Data to plot.
    obs_param : str, optional
        obs parameter (zp/mean_wave). The default is 'zp'.
    unit_obs_param : str, optional
        unit for the obs param. The default is 'mmag'.
    vary : str, optional
        y-axis variable. The default is 'sigma_atmos_param'.
    err_rel : bool, optional
        To plot relative errors or not. The default is False.
    atmos_params : list(str), optional
        List of atmos params. The default is ['airmass','ozone','aerosol','pwv'].
    limy : list(list(float)), optional
        y-axis limits for the plot. The default is [[0.,0.05],[0.,100.],[0.,0.0055],[0.,0.5]].
    auxtel_data : list(float), optional
        Auxtel typical values. The default is [3.e-3,20,5.e-3,0.2].
    plotDir : str, optional
        Output dir for the plots. The default is ''.

    Returns
    -------
    None.

    """
    from sn_analysis.sn_atmos_tools import add_index_band,add_legend
    
    res = add_index_band(res)
    if err_rel:
        res['sigma_atmos_param']/=res['atmos_param_value']/100.
        #limy = [0.,30.]*4
    
    auxtel_mes = dict(zip(atmos_params,auxtel_data))
    bands_atm = dict(zip(atmos_params,
                 ['grizy','gri','grizy','izy']))
    lstyle = dict(zip([1.2,2.0],['solid','dotted']))
    
    limy = dict(zip(atmos_params,limy))
    unit = dict(zip(atmos_params,['','[DU]','','[mm]']))

    sigmas = res['sigma_obs_param'].unique()
    
    markers = ['o','s','P','h']
    colors = ['m','r','b','g']
    
    mm = dict(zip(sigmas,markers))
    ccolors = dict(zip(sigmas,colors))
    
    for atm in atmos_params:
        idx = res['atmos_param'] == atm
        sela = res[idx]
        fig, ax = plt.subplots(figsize=(12,8))
        fig.subplots_adjust(right=0.85)
        airmass = sela['airmass'].unique()
        
        for airm in airmass:
            idx = sela['airmass'] == airm
            idx &= sela['band'].isin(list(bands_atm[atm]))
            selb = sela[idx]
            selb = selb.sort_values(by=['index','sigma_obs_param'])
            sigmas = selb['sigma_obs_param'].unique()
            
            for sig in sigmas:
                thelab = None
                if  airm == 1.2:
                    thelab = '$\sigma_{'+obs_param+'}$='+'{}'.format(sig)+' '+unit_obs_param
                
                idx = selb['sigma_obs_param'] == sig
                selc = selb[idx]
                ax.plot(selc['band'],selc[vary],
                        color=ccolors[sig],linestyle=lstyle[airm],
                        marker=mm[sig],mfc='None',label=thelab)
            
        ax.grid(visible=True)
        ax.set_ylim(limy[atm])
    
        if not err_rel:
            sig_atm = '$\sigma_{'+atm+'}$ '+format(unit[atm])
        else:
            sig_atm = '$\\frac{\sigma_{'+atm+'}}{<'+atm+'>}$ [%]'
            
        ax.set_ylabel(sig_atm)
        ax.set_xlabel('band')
        ax.legend(loc='upper left',
              bbox_to_anchor=(-0.1, 1.15), ncol=4, frameon=False, fontsize=15)
       
        # add airmass legend
        add_legend(ax, airmass)
        
        # add auxtel typical measurements
        
        xmin,xmax = ax.get_xlim()
        vmes = auxtel_mes[atm]
        ax.plot([xmin,xmax],[vmes]*2,linestyle='dashed',color='k')
        ax.set_xlim([xmin,xmax])
        vunit = unit[atm].split('[')[-1].split(']')[0]
        if not err_rel:
            ttext = '$\sigma_{'+atm+'}$='+'{}'.format(vmes)+ ' {}'.format(vunit)
        else:
            import numpy as np
            ttext = '$\\frac{\sigma_{'+atm+'}}{<'+atm+'>}$='+'{}'.format(np.round(vmes,1))+ ' %'
        ax.text(1.02*xmax,vmes,ttext,fontsize=12)
    
        if plotDir != '':
            outName = '{}/sigma_{}_{}.png'.format(plotDir,atm,int(err_rel))
            plt.savefig(outName)