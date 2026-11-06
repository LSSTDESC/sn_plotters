#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 09:33:06 2025

@author: philippe.gris@clermont.in2p3.fr
"""
from . import plt
import pandas as pd
import numpy as np
from datetime import datetime
from astropy.time import Time
from sn_tools.sn_obs import season
from sn_tools.sn_io import checkDir
import healpy as hp
import matplotlib
from sn_tools.sn_utils import multiproc
import os
import numpy.lib.recfunctions as rf

plt.rcParams["axes.labelsize"] = "large"
plt.rcParams["axes.linewidth"] = 2.0
plt.rcParams["xtick.major.size"] = 8
plt.rcParams["ytick.major.size"] = 8
plt.rcParams["ytick.minor.size"] = 5
plt.rcParams["xtick.labelsize"] = "large"
plt.rcParams["ytick.labelsize"] = "large"

plt.rcParams["figure.figsize"] = (12, 8)
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.titleweight'] = 'bold'
# plt.rcParams['axes.facecolor'] = 'blue'
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'
# the line width around the marker symbol
plt.rcParams['lines.markeredgewidth'] = 0.3
plt.rcParams['lines.markersize'] = 5  # markersize, in points
plt.rcParams['grid.alpha'] = 0.75  # transparency, between 0.0 and 1.0
plt.rcParams['grid.linestyle'] = '-'  # simple line
plt.rcParams['grid.linewidth'] = 0.4  # in points
plt.rcParams['font.size'] = 13


def add_info(df):
    """
    Function to add infos about the df

    Parameters
    ----------
    df : pandas df
        Data to process.

    Returns
    -------
    df : pandas df
        Output data.

    """

    # add the date
    dfb = pd.DataFrame({'year': df['year'].to_list(),
                        'month': df['month'].to_list(),
                        'day': df['day'].to_list()})

    df['date'] = pd.to_datetime(dfb).to_list()

    # estimate Nvisits corresponding to obs_time
    df['nvisits'] = (df['obs_duration [h]']*3600/30).astype(int)

    # add seasons
    obs = df.to_records(index=False)

    targets = np.unique(obs['target'])

    res = None
    for tt in targets:
        idx = obs['target'] == tt
        sel = obs[idx]
        sel_seas = season(sel, season_gap=20, mjdCol='mjd')
        if res is None:
            res = sel_seas
        else:
            res = np.concatenate((res, sel_seas))

    df = pd.DataFrame.from_records(res)
    df = df.rename(columns={'period': 'season'})

    mjd_min = df['mjd'].min()
    df['night'] = (df['mjd']-mjd_min+1).astype(int)

    return df


def plot(df, fields=['COSMOS']):
    """
    Function to plot observing time vs time

    Parameters
    ----------
    df : pandas df
        Data to process.
    fields: list(str)
        List of fields to display.

    Returns
    -------
    None.

    """

    mjd_min = df['mjd'].min()
    ymax = df['nvisits'].max()+10

    for tt in fields:
        fig, ax = plt.subplots(figsize=(14, 8))
        fig.suptitle(tt)
        idx = df['target'] == tt
        sel = pd.DataFrame(df[idx])
        """
        ax.plot(sel['date'], sel['obs_duration [h]'],
                marker='.', linestyle='None')
        """
        ax.plot(sel['date'], sel['nvisits'],
                marker='.', linestyle='None', color='b')

        seasons = sel['season'].unique()

        for seas in seasons:
            ii = sel['season'] == seas
            selb = sel[ii]
            ax.plot(selb['date'], selb['nvisits'],
                    marker='.', linestyle='None')

    for nn in range(11):
        xv = mjd_min+(nn*365)
        ti = Time('{}'.format(xv), format='mjd')
        print('bbb', ti.ymdhms)
        ddate = datetime(ti.ymdhms[0], ti.ymdhms[1], ti.ymdhms[2])
        ax.plot([ddate]*2, [0, ymax], color='r', lw=2)

    ymin = 0.

    ax.set_ylim([ymin, ymax])
    ax.grid(visible=True)
    ax.set_xlabel(r'Time [year]')
    ax.set_ylabel(r'$N_{visits}$')


def ana_season(sel_df):
    """
    Function to analyze a season

    Parameters
    ----------
    sel_df : pandas df
        Data to process.

    Returns
    -------
    pandas df
        output data.

    """
    tti = np.arange(0., 4.5, 0.02)

    var = 'obs_duration [h]'
    r = []
    for i in range(len(tti)-1):
        io = sel_df[var] >= tti[i]
        io &= sel_df[var] < tti[i+1]
        ssel = sel_df[io]
        sl = ssel['mjd'].max()-ssel['mjd'].min()
        r.append((np.mean([tti[i], tti[i+1]]), sl))

    res = np.rec.fromrecords(r, names=['observable_time [h]', 'season_length'])

    return pd.DataFrame.from_records(res)


def smooth_It(vals, xvar, yvar, kk=3):
    """
    Function to smooth a set of data

    Parameters
    ----------
    vals : array
        Data to process.
    xvar : str
        x-axis variable.
    yvar : str
        y-axis variable.
    kk : int, optional
        smoothing parameter. The default is 3.

    Returns
    -------
    xnew : array
        xnew.
    spl_smooth : array
        y smoothed values.

    """

    from scipy.interpolate import make_interp_spline, UnivariateSpline
    xmin, xmax = np.min(vals[xvar]), np.max(vals[xvar])
    xnew = np.linspace(xmin, xmax, 100)
    spl = make_interp_spline(
        vals[xvar], vals[yvar], k=kk)  # type: BSpline
    spl = UnivariateSpline(vals[xvar], vals[yvar], k=kk)
    spl.set_smoothing_factor(0.5)
    spl_smooth = spl(xnew)
    return xnew, spl_smooth


def fa(x):
    return x*3600./30.


def fb(x):

    return x*3600./30.


def plot_season_length(res_season):
    """
    Function to plot season length vs obs time/Nvisits

    Parameters
    ----------
    res_season : pandas df
        Data to plot.

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots()
    # axb = ax.twiny()
    dfb = pd.DataFrame()
    tta = ['COSMOS', 'XMM_LSS', 'ELAISS1', 'ECDFS', 'EDFS_a']
    lcol = ['g', 'm', 'r', 'orange', 'b']
    lst = ['solid', 'dotted', 'dashdot', 'dashed', 'dashdot']
    labs = ['COSMOS', 'XMM-LSS', 'ELAIS', 'CDFS', 'EDFS$_{a,b}$']

    dcol = dict(zip(tta, lcol))
    dlst = dict(zip(tta, lst))
    dlabs = dict(zip(tta, labs))

    for tt in tta:
        idx = res_season['target'] == tt
        idx &= res_season['season_length'] > 50
        sel_seas = res_season[idx]
        # sel_seas['nvisits'] = (sel_seas['observable_time [h]']*3600/30).astype(int)
        dfb = pd.concat((dfb, sel_seas))

        ax.plot(sel_seas['observable_time [h]'],
                sel_seas['season_length'], marker='None',
                label=dlabs[tt], color=dcol[tt], linestyle=dlst[tt], lw=2)
        # axb.plot(sel_seas['nvisits'],
        #         sel_seas['season_length'])
        # x, y = smooth_It(sel_seas, 'observable_time [h]', 'season_length')
        # ax.plot(x, y, marker='.', label=tt)

    ax.grid(visible=True)

    ymin = dfb['season_length'].min()-5.
    ymax = dfb['season_length'].max()+5.
    ax.plot([0, 4], [180.]*2, color='k', lw=2)
    ax.set_xlim([0, 4])
    ax.set_ylim([ymin, ymax])
    ax.legend()
    ax.set_xlabel(r'Observable time [h]')
    ax.set_ylabel(r'Max season length [days]')
    ax.text(0.5, 185., 'season length = 6 months')
    secax = ax.secondary_xaxis('top', functions=(fa, fb))
    secax.set_xlabel('N$_{visits}$')

class LSST_plus_AuxTel_visu:
    def __init__(self,obs,dbName,nside=32,figDir='None',
                 movieDir='None',ffmpeg_cmd='ffmpeg',year=1):
        """
        class to visualize LSST and AuxTel surveys

        Parameters
        ----------
        obs : array
            observations.

        Returns
        -------
        None.

        """
        
        idx = obs['year'] == year
        self.obs = obs[idx]
        self.year=year
        self.npixels = hp.nside2npix(nside)
        ra_min ,ra_max = np.min(obs['RA']),np.max(obs['RA'])
        dec_min ,dec_max = np.min(obs['Dec']),np.max(obs['Dec'])
        
        self.coord_lim = {}
        
        self.coord_lim['RA'] = [ra_min,ra_max]
        self.coord_lim['Dec'] = [dec_min,dec_max]
        
        # prepare output directory
        self.plotDir = 'None'
        if figDir != 'None':
            self.plotDir = '{}/{}/year_{}'.format(figDir,dbName,year)
            checkDir(self.plotDir)
        
        self.dbName = dbName
        
        self.movieDir = movieDir
        if self.movieDir != 'None':
            self.ffmpeg_cmd = ffmpeg_cmd
            checkDir(self.movieDir)
            
    def __call__(self,nights=[2],nproc=8):
           
        params = {}
        
        if nproc > 1:
            from sn_tools.sn_utils import multiproc
            multiproc(nights,params,self.process,nproc)
        else:
            self.process(nights,params)
           
    def process(self,nights, params, j=0, output_q=None):
        
        for nnight in nights:
            print('processing night',nnight)
            night = nnight+365*(self.year-1)
            idx = self.obs['night']==night
            obs_night = self.obs[idx]
            if len(obs_night) == 0:
                continue
            idxb = self.obs['night'] < night
            obs_bef = self.obs[idxb]
            
            obs_bef = rf.append_fields(obs_bef,'count',
                                       [1]*len(obs_bef),
                                       dtypes='<i8')
            self.loop_obs_night(obs_night,obs_bef,night)
        
            del obs_bef
            del obs_night
            
            if self.movieDir != 'None':
                self.makeMovie(night)
        
        if output_q is not None:
            return output_q.put({j: 1})
        else:
            return 1
           
           
           
    def loop_obs_night(self,obs,obs_bef,night):
        """
        Method to loop on observations
 
        Returns
        -------
        None.
 
        """
        
        colNames = ['MJD','dist_obs_calib',
                    'deltaT [min]','slew_time [s]']
        
        df_pl = pd.DataFrame()
        for i in range(1,len(obs)):
            r = []
            mjd = obs[i]['mjd']
            if i > 0:
                dd_cur = self.get_vals(obs[i])
                dd_bef = self.get_vals(obs[i-1])
               
                deltaT = 24*60*(dd_cur['delta_mjd'])
                slew_time = dd_bef['slewTime']
                dist = dd_cur['delta_r']
                r.append((mjd,dist,deltaT,slew_time))
        
            else:
                r.append((mjd,0,0,0))
                
            df_plb = pd.DataFrame(r,columns=colNames)
            df_pl = pd.concat((df_pl,df_plb))
        
            obs_pl = obs[:i+1]
            
            rb = [2]*len(obs_pl)
            rb[-1] = 3
            if len(obs_pl) == 1:
                rb = [3]
        
            obs_pl = rf.append_fields(obs_pl,'count',rb)
            
            self.plot_night(obs_bef,obs_pl,df_pl,night,i)
    
        del df_pl
           
    def get_vals(self,obs):
        
        cols = ['RA','Dec','mjd','slewTime',
                'RA_deg','Dec_deg','delta_mjd','delta_r']
       
        dd = dict(zip(cols,[obs[col] for col in cols]))
       
        return dd
           
       
    def plot_night(self,obs_bef,obs,df_pl,night,inum):
        
        fig = plt.figure(figsize=(15, 12))
        fig.subplots_adjust(top=0.99)
        nnight = night-365*(self.year-1)
        fig.suptitle('year {} - night {}'.format(self.year,nnight))
        ax_a = fig.add_axes([0.1, 0.50, 0.34, 0.40])
        ax_b = fig.add_axes([0.1, 0.01, 0.3, 0.35])
        ax_c = fig.add_axes([0.5, 0.65, 0.35, 0.25])
        ax_d = fig.add_axes([0.5, 0.4, 0.35, 0.25])
        ax_e = fig.add_axes([0.5, 0.15, 0.35, 0.25])
        
        self.plot_ra_dec(ax_a,obs,obs_bef)
        self.plot_Mollview(ax_b,obs[['healpixID','count']],
                           obs_bef[['healpixID','count']],
                           fig,night)
        """
        self.plot_info(ax_c,df_pl,
                       varx='MJD',legx='MJD',
                       vary='dist_obs',legy='dist [deg]',label='obs')
        self.plot_info(ax_c,df_pl,
                       varx='MJD',legx='MJD',
                       vary='dist_calib',
                       legy='dist [deg]',color='b',label='calib star',
                       lstyle='dashed')
        """
        self.plot_info(ax_c,df_pl,
                      varx='MJD',legx='MJD',
                      vary='dist_obs_calib',legy='dist [deg]',
                      color='r',label='obs/calib star',lstyle='solid')
        
        self.plot_info(ax_d,df_pl,
                       varx='MJD',legx='MJD',
                       vary='deltaT [min]',legy='$\Delta t$ [min]',
                       label='$\Delta t$')
        
        self.plot_info(ax_e,df_pl,
                       varx='MJD',legx='MJD',
                       vary='slew_time [s]',legy='LSST slew time [s]',
                       color='m',lstyle='dashed',
                       label='')
        """
        self.plot_info(ax_e,df_pl,
                       varx='MJD',legx='MJD',
                       vary='slewing_speed [deg/s]',
                       legy='slewing \n speed [deg/s]')
        """
        ax_c.set_xticklabels([])
        ax_d.set_xticklabels([])
        
        ax_a.legend(bbox_to_anchor=(0.85, 1.15), ncol=2,
             fontsize=10, frameon=False)
        ax_c.legend(bbox_to_anchor=(0.85, 1.15), ncol=3,
             fontsize=10, frameon=False)
        ax_d.legend(bbox_to_anchor=(0.95, 0.97), ncol=2,
             fontsize=10, frameon=False)
        
        
        if self.plotDir != 'None':
            outDir = '{}/night_{}'.format(self.plotDir,nnight)
            checkDir(outDir)
            fName = '{}/fig_{}.jpg'.format(outDir,str(inum).zfill(3))
            plt.savefig(fName)
            plt.close(fig)
        else:
            plt.show()
        
    def plot_ra_dec(self,ax,obs,obs_bef,varx='RA',legx='RA [deg]',
                    vary='Dec',legy='Dec [deg]'):
        
        
        plt.sca(ax)
        msize=5
        if len(obs_bef) > 0:
            ax.plot(obs_bef[varx],obs_bef[vary],
                color='grey',marker='o',
                markersize=msize,linestyle='None')
            ax.plot(obs_bef['RA_deg'],obs_bef['Dec_deg'],
                color='lightskyblue',marker='*',markersize=8,
                linestyle='None')
       
        
        n = len(obs)-1
        ax.plot(obs[:n][varx],obs[:n][vary],
                color='k',marker='o',markersize=msize,
                linestyle='None',label='LSST pointings')
  
        ax.plot(obs[n][varx],obs[n][vary],
                color='r',marker='o',
                markersize=msize,linestyle='None')
        
        ax.plot(obs[n]['RA_deg'],obs[n]['Dec_deg'],
                color='b',marker='*',markersize=8,
                linestyle='None',label='calib star')
        
        ax.set_xlim([self.coord_lim[varx][0],self.coord_lim[varx][1]])
        ax.set_ylim([self.coord_lim[vary][0],self.coord_lim[vary][1]])
        
        ax.set_xlabel(r'{}'.format(legx))
        ax.set_ylabel(r'{}'.format(legy))
        ax.grid(visible=True)
        
    def plot_Mollview(self, axa,obs,obs_bef,fig, night,var='count'):
        """
        Method to display a Mollweid view
 
        Parameters
        --------------
        pixels: pandas df
          data to plots
        axa: matplotlib axis
          axis to use for the plot
        fig: matplotlib figure
          figure to use for plot
        night: int
          night number
        """
        
        
        plt.sca(axa)
 
        xmin = 0.99999
        xmax = np.max([np.max(obs[var]), 1])
 
        norm = plt.cm.colors.Normalize(xmin, xmax)
        
        cmap = matplotlib.colors.ListedColormap(["grey", "black", "red"])
        cmap.set_under('w')
        """
        cmap = plt.get_cmap('jet', xmax)
        cmap.set_under('w')
        """
    
        hpxmap = np.zeros(self.npixels, dtype=int)
        hpxmap = np.full(hpxmap.shape, 0)
        hpxmap[obs_bef['healpixID']] = obs_bef[var]
        hpxmap[obs['healpixID']] = obs[var].astype(int)
 
    
        # print('hello ',xmin,xmax)
 
        hp.mollview(hpxmap, nest=True, cmap=cmap,
                    min=xmin, max=xmax, norm=norm, cbar=False,
                    title='', hold=True, badcolor='white', xsize=1600)
 
        hp.graticule()
 
    def plot_info(self,ax,df_pl,
                      varx='MJD',legx='MJD',
                      vary='dist',legy='dist [deg]',color='k',
                      label=None,lstyle='solid'):
        
        ax.plot(df_pl[varx],df_pl[vary],color=color,
                label=label,linestyle=lstyle)
        
        ax.set_ylabel(r'{}'.format(legy),rotation=270,labelpad=50)
        
        ax.yaxis.tick_right()
        ax.grid(visible=True)
        ax.yaxis.set_label_position("right")
        ax.set_xlabel(r'{}'.format(legx))
        
        
    def makeMovie(self,night):
         """
         Method to make a movie from jpg files
         """
         nnight = night-365*(self.year-1)
         dirFigs = '{}/night_{}'.format(self.plotDir,nnight)
         movieDir = '{}/year_{}'.format(self.movieDir,self.year)
         checkDir(movieDir)
         # ffmpeg_cmd = '../../ffmpeg/ffmpeg-4.3.1-i686-static/ffmpeg'
         
         cmd = '{} -r 2 -f image2 -s 1920x1080 -i {}/fig_%03d.jpg -vcodec libx264 -crf 25  -pix_fmt yuv420p {}/{}_{}.mp4 -y'.format(
             self.ffmpeg_cmd, dirFigs, movieDir, self.dbName,nnight)
         os.system(cmd)