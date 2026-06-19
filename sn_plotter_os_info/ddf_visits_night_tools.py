#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 10:01:53 2025

@author: philippe.gris@clermont.in2p3.fr
"""
import numpy as np
from . import plt, filtercolors, filtermarkers
import pandas as pd
from sn_analysis.sn_calc_plot import bin_it
import operator as op
from sn_plotter_analysis.sn_analyser_tools import clean_level


def analyze_simu_exp(data):
    """
    Function to analyze sim+expected data

    Parameters
    ----------
    data : pandas df
        Data to process.

    Returns
    -------
    None.

    """

    # add new columns as diff

    for vv in ['nvisits', 'u', 'g', 'r', 'i', 'z', 'y']:
        data['diff_{}'.format(vv)] = data['{}_exp'.format(vv)] - data[vv]
        data['ratio_{}'.format(vv)] = data['{}_exp'.format(vv)] / data[vv]

    print(data['ratio_u'])
    data = data.replace([np.inf, -np.inf], 0)
    print(data['ratio_u'])

    # plot stat results
    ana_plot_stat(data)

    # analysis of cases when thee number of visits exceeds expectation
    plot_stat_visits_vs_exp(data, op.lt, '<')

    # plot obs time vs night
    plot_obs_time_night(data)

    plt.show()


def plot_obs_time_night(data):
    """
    Function to plot obs time and nddf per night of observation

    Parameters
    ----------
    data : pandas df
        Data to process.

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots(nrows=2, figsize=(12, 8))
    fig.subplots_adjust(hspace=0)
    print(data.columns)
    fig.suptitle(data['dbName'].unique()[0])
    obs_time = data.groupby(['night']).apply(
        lambda x: pd.DataFrame(
            {'obs_time [h]': [x['nvisits'].sum()*30./3600.],
             'nddf': [len(x['field'].unique())]
             })).reset_index()

    ax[0].plot(obs_time['night'], obs_time['nddf'], 'k.', ms=5)

    ax[0].set_xticklabels([])
    ax[0].set_ylabel(r'N$_{DDF}$')

    ax[1].plot(obs_time['night'], obs_time['obs_time [h]'], 'k.', ms=5)
    print(obs_time)

    ax[1].set_xlabel(r'night')
    ax[1].set_ylabel(r'DDF obs. time [h]')

    night_max = obs_time['night'].max()
    tdays = np.arange(1, night_max, 365.)

    colors = ['blue', 'orange', 'green', 'crimson', 'lightsalmon',
              'grey', 'yellow', 'magenta', 'purple', 'cyan']
    print(len(tdays))

    for i in range(2):
        ax[i].grid(visible=True)
        ax[i].set_xlim([0, night_max])
        ax[i].set_ylim([0, None])
        # ax[i].axvspan(1, 365, facecolor='red', alpha=0.25)
        for j in range(len(tdays)):
            ax[i].axvspan(tdays[j], tdays[j]+365.,
                          facecolor=colors[j], alpha=0.25)


def plot_stat_visits_vs_exp(data, ope, opevalue, selval=0.,
                            field='DD:COSMOS', bands='grizy',
                            bins=np.arange(0, 1.1, 0.1)):
    """
    Function to plots diff exp/obs per night

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    idx = ope(data['diff_nvisits'], selval)
    idx &= data['field'] == field
    sel = data[idx]

    dd = {}
    for b in bands:
        vvar = 'ratio_{}'.format(b)
        rb = bin_it(sel, vvar, bins=bins, norm_factor=1, outvar='frac')
        rb['frac'] /= rb['frac'].sum()
        dd[b] = rb

    fig, ax = plt.subplots(figsize=(12, 8))
    fig.subplots_adjust(hspace=0)
    dbName = sel['dbName'].unique()[0]
    year = sel['year'].unique()[0]
    field = sel['field'].unique()[0]
    figtit = '{} - {} \n year {}'.format(dbName, field, year)

    figtit += '- $\\frac{N_{visits}^{obs}}{N_{visits}^{exp}}$'+opevalue+'1'

    fig.suptitle(figtit)

    for key, vals in dd.items():
        ax.plot(vals['ratio_{}'.format(key)], 100.*vals['frac'],
                label=key,
                color=filtercolors[key],
                marker=filtermarkers[key], mfc='None')

    ax.grid(visible=True)
    # ax.set_ylabel(r'Fraction of nights [%]')
    laby = r'N$_{nights}$ [%]'
    ax.set_ylabel(laby)
    ax.set_xlabel(r'$\frac{N_{visits}^{obs}}{N_{visits}^{exp}}$')
    xmin = np.min(bins)
    ax.set_xlim([xmin, None])
    ax.set_ylim([0, None])

    ax.legend()
    fig.tight_layout()
    """
    ax.hist(sel['ratio_nvisits'], histtype='step', bins=20)

    for b in 'ugrizy':
        fig, ax = plt.subplots()
        fig.suptitle('{}-band'.format(b))
        vvar = 'ratio_{}'.format(b)
        ax.hist(sel[vvar], histtype='step', bins=20)

    plt.show()
    """


def get_nights(data, colName, selval, bands='grizy'):
    """
    getting nights corresponding to specific filter alloc configs

    Parameters
    ----------
    data : pandas df
        Data to process.
    colName : str
        colName to use.
    selval : float
        selection values.
    bands : str, optional
        bands to consider. The default is 'grizy'.

    Returns
    -------
    None.

    """

    idx = True

    for b in bands:
        idx &= np.abs(data['{}{}'.format(colName, b)]-selval) < 0.01

    sel = data[idx]

    nights = np.unique(sel['night'])

    return nights


def ana_plot_stat(data):
    """
    Function to make some stat and to make plots

    Parameters
    ----------
    data : pandas df
        Data to process.

    Returns
    -------
    nights: list(int)
      list of corresponding nights

    """

    res_stat = data.groupby(['field', 'season', 'DD_type', 'dbName']).apply(
        lambda x: stat_simu_exp(x)).reset_index()

    print(res_stat)

    plot_stat(res_stat)


def stat_simu_exp(grp):
    """
    Function to get some stat on simu+expected os data

    Parameters
    ----------
    grp : pandas df
        Data to process.

    Returns
    -------
    rr : pandas df
        output data.

    """

    dict_frac = {}
    nnights = len(grp)
    for vv in ['nvisits', 'u', 'g', 'r', 'i', 'z', 'y']:
        myvar = 'diff_{}'.format(vv)
        idxa = grp[myvar] >= 1
        idxb = grp[myvar] <= -1
        idxc = np.abs(grp[myvar]) < 0.5
        dict_frac['{}_missing'.format(vv)] = [len(grp[idxa])/nnights]
        dict_frac['{}_excess'.format(vv)] = [len(grp[idxb])/nnights]
        dict_frac['{}_perfect'.format(vv)] = [len(grp[idxc])/nnights]

    rr = pd.DataFrame.from_dict(dict_frac)

    return rr


def plot_stat(data, prefix='DD:'):
    """
    Function to plot stat results

    Parameters
    ----------
    data : pandas df
        Data to process.
    prefix : str, optional
        prefix for DD fieldnames. The default is 'DD:'.

    Returns
    -------
    None.

    """
    print(data.columns)
    fields = data['field'].unique()

    categ = ['perfect', 'missing', 'excess']
    value = ['=', '<', '>']
    fields = ['COSMOS', 'XMM_LSS', 'ELAISS1', 'ECDFS', 'EDFS_a', 'EDFS_b']
    fields = list(map(lambda x: prefix + x, fields))
    colors = ['r', 'b', 'k', 'orange', 'm', 'g']
    marks = ['o', 's', '*', '^', 'v', '>']
    dict_col = dict(zip(fields, colors))
    dict_mark = dict(zip(fields, marks))
    dict_ls = dict(zip(['UD', 'DF'], ['solid', 'dashed']))

    for i, vval in enumerate(categ):
        fig, ax = plt.subplots(figsize=(14, 8))
        figtit = data['dbName'].unique()[0]
        figtit += r' - $\frac{N_{visits}^{exp}}{N_{visits}^{simu}}$' + \
            value[i]+'1'
        fig.suptitle(figtit)
        fig.subplots_adjust(right=0.85)
        for field in fields:
            idx = data['field'] == field
            sel = data[idx]
            dd_type = sel['DD_type'].unique()[0]
            ax.plot(sel['season'], sel['nvisits_{}'.format(vval)],
                    color=dict_col[field], marker=dict_mark[field],
                    linestyle=dict_ls[dd_type],
                    mfc='None', ms=10, label=field.split(prefix)[-1])

        ax.grid(visible=True)
        ax.set_xlabel(r'season')
        """
        ax.set_ylabel(
            r'$\frac{N_{visits}^{exp}}{N_{visits}^{simu}}$'+value[i]+'1')
        
        """
        laby = r'N$_{nights}$ [%]'
        ax.set_ylabel(laby)
        ax.legend(bbox_to_anchor=(0.99, 0.7),
                  ncol=1, fontsize=15, frameon=False)


def ana_seq_multi(toproc, params, j=0, output_q=None):
    """
    Analysis function using multiprocessing

    Parameters
    ----------
    toproc : list(str)
        List of OS to process.
    params : dict
        parameters.
    j : int, optional
        internal tag for multiprocessing. The default is 0.
    output_q : multiprocessing queue, optional
        Where to put the data. The default is None.

    Returns
    -------
    pandas df
        Analyzed data.

    """

    timescale = params['timescale']
    df = params['data']

    idx = df['dbName'].isin(toproc)

    sel = pd.DataFrame(df[idx])

    del df

    res = ana_seq(sel, timescale)

    if output_q is not None:
        return output_q.put({j: res})
    else:
        return res


def ana_seq(df, timescale='year'):
    """
    Function to analyze DDF sequences

    Parameters
    ----------
    df : pandas df
        Data to process.
    timescale : str, optional
        Timescale. The default is 'year'.

    Returns
    -------
    dfd : pandas df
        Output data.

    """

    ccols_m = ['field', timescale, 'dbName']
    ccols = ccols_m+['seq_tot']

    dfb = df.groupby(ccols)[ccols].apply(
        lambda x: get_nvisits(x)).reset_index()
    dfb = clean_level(dfb)

    bands = 'ugrizy'
    colsb = ccols+list(bands)
    # dfb = dfb.merge(df[colsb], left_on=ccols,right_on=ccols,suffixes=['',''])

    for b in bands:
        ccols = ccols_m+[b]+['seq_tot']
        ccob = 'nnights_{}'.format(b)
        dfe = df.groupby(ccols)[ccols].apply(
            lambda x: get_nvisits_band(x, b, ccob)).reset_index()

        dfe = clean_level(dfe)

        dfb = dfb.merge(dfe, left_on=ccols_m+['seq_tot'],
                        right_on=ccols_m+['seq_tot'], suffixes=['', ''])

    ccols = ccols_m+['night']

    dfc = df.groupby(ccols_m)[ccols].apply(
        lambda x: get_nnights(x)).reset_index()
    dfc = clean_level(dfc)
    dfd = dfb.merge(dfc, left_on=ccols_m, right_on=ccols_m, suffixes=['', ''])

    dfd['seq_frac'] = 100.*dfd['nnights']/dfd['nnights_year']

    return dfd


def get_nvisits(grp, thevar='nnights'):
    """
    Function to estimate the number of nights corresponding to a DDF sequence

    Parameters
    ----------
    grp : pandas df
        Data to process.
    thevar : str, optional
        output col name. The default is 'nnights'.

    Returns
    -------
    res : pandas df
        output data.

    """

    dd = {}

    dd[thevar] = [len(grp)]

    res = pd.DataFrame.from_dict(dd)

    res[thevar] = res[thevar].astype(int)
    return res


def get_nvisits_band(grp, thevar, thevar_name='nnights'):
    """
    Function to get the number of visits per band


    Parameters
    ----------
    grp : pandas df
        Data to process.
    thevar : str
        col to process.
    thevar_name : str, optional
        col output name. The default is 'nnights'.

    Returns
    -------
    res : pandas df
        output result.

    """

    dd = {}

    idx = grp[thevar] > 0
    sel = grp[idx]

    dd[thevar_name] = [len(sel)]

    res = pd.DataFrame.from_dict(dd)

    res[thevar_name] = res[thevar_name].astype(int)
    return res


def get_nnights(grp, thevar='nnights_year'):
    """
    Function to estimate the total number of nights

    Parameters
    ----------
    grp : pandas df
        Data to process.
    thevar : str, optional
        output col name. The default is 'nnights_year'.

    Returns
    -------
    res : pandas df
        Result.

    """

    dd = {}
    nights = grp['night'].unique()
    dd[thevar] = [len(nights)]

    res = pd.DataFrame.from_dict(dd)

    res[thevar] = res[thevar].astype(int)

    return res


def calc_summary(grp, col='y'):
    """
    Function to extract some result

    Parameters
    ----------
    grp : pandas df
        Data to process.
    col : str, optional
        col name to select. The default is 'y'.

    Returns
    -------
    rr : pandas df
        output result.

    """

    idx = grp[col] > 0
    sel = grp[idx]
    selb = sel.sort_values(by=['seq_frac'], ascending=False)

    rr = selb[['seq_tot', 'seq_frac', col]][:1]
    rr['nvisits_{}'.format(col)] = selb[col].sum()
    nnights_band = selb['nnights_{}'.format(col)].sum()
    rr['frac_{}'.format(col)] = sel['seq_frac'].sum()
    # correct to get the fraction of seq corresponding to the band
    nnights_year = selb['nnights_year'][:1]
    rr['seq_frac'] *= nnights_year/nnights_band

    rr = rr.rename(columns={'seq_tot': 'seq_tot_{}'.format(col),
                            'seq_frac': 'seq_frac_{}'.format(col)})
    return rr


def summary_seq(grp):
    """
    function to estimate summary results

    Parameters
    ----------
    grp : pandas df
        Data to process.

    Returns
    -------
    res : pandas df
        output data.

    """

    rr = calc_summary(grp, 'y')

    bb = calc_summary(grp, 'u')

    res = rr.merge(bb, how='cross')

    return res


def get_ratios(grp, df_orig, band='y'):
    """
    Function to extract the ratios of the number of visits per band wrt ref

    Parameters
    ----------
    grp : pandas df
        Data to process.
    df_orig : pandas df
        Data to extract info from.
    band : str, optional
        band considered. The default is 'y'.

    Returns
    -------
    sel_test : pandas df
        The result.

    """

    dbName = grp.name[0]
    target_name = grp.name[1]
    year = grp.name[2]

    idx = df_orig['dbName'] == dbName
    idx &= df_orig['field'] == target_name
    idx &= df_orig['year'] == year
    idx &= df_orig[band] > 0

    sel_orig = pd.DataFrame(df_orig[idx])

    del df_orig

    nnights_y = len(sel_orig)

    b_ref = get_ref_sequence(grp, band)

    b_ref = clean_level(b_ref)
    sel_test = sel_orig.merge(b_ref, how='cross')

    sel_test['diff_nvisits'] = sel_test['nvisits']-sel_test['nvisits_ref']
    bands = 'grizy'
    for b in bands:
        sel_test['ratio_{}'.format(
            b)] = sel_test[b]/sel_test['{}_ref'.format(b)]

    sel_test = clean_level(sel_test)

    return sel_test


def get_ref_sequence(grp, band):
    """
    Function to extract the reference sequence as a df

    Parameters
    ----------
    grp : pandas df
        Data to process.
    band : str
        band to consider for the ref sequence.

    Returns
    -------
    b_ref : pandas df
        Reference df sequence.

    """

    seq = grp['seq_tot_{}'.format(band)].values[0]

    spl = seq.split('-')

    bands = [sp[-1] for sp in spl]
    bands_ref = list(map(lambda el: el+'_ref', bands))
    nv_ref = list(map(int, [sp[:-1] for sp in spl]))
    nv_ref = list(map(lambda el: [el], nv_ref))
    dd = dict(zip(bands_ref, nv_ref))
    b_ref = pd.DataFrame.from_dict(dd)

    b_ref['nvisits_ref'] = b_ref[bands_ref].sum(axis=1)

    return b_ref


def get_stat_indiv(grp, df_orig, band='y'):
    """
    Function to grab the number of nights corresponding to a sequence 

    Parameters
    ----------
    grp : pandas df
        Data to process.
    df_orig : pandas df
        Data to process.
    band : str, optional
        filter for the sequence. The default is 'y'.

    Returns
    -------
    res : pandas df
        Result.

    """

    dbName = grp.name[0]
    target_name = grp.name[1]
    year = grp.name[2]

    idx = df_orig['dbName'] == dbName
    idx &= df_orig['field'] == target_name
    idx &= df_orig['year'] == year
    idx &= df_orig[band] > 0

    sel_orig = pd.DataFrame(df_orig[idx])

    del df_orig

    b_ref = get_ref_sequence(grp, band)

    b_ref = clean_level(b_ref)
    sel_test = sel_orig.merge(b_ref, how='cross')

    sel_test['diff_nvisits'] = sel_test['nvisits']-sel_test['nvisits_ref']
    bands = 'grizy'
    for b in bands:
        sel_test['ratio_{}'.format(
            b)] = sel_test[b]/sel_test['{}_ref'.format(b)]

    # three types of nights
    dd = {}
    dd['frac_equal'] = [100.*get_val(sel_test, 'diff_nvisits', op.eq, 0)]
    dd['frac_plus'] = [100.*get_val(sel_test, 'diff_nvisits', op.gt, 0)]
    dd['frac_minus'] = [100.*get_val(sel_test, 'diff_nvisits', op.lt, 0)]

    res = pd.DataFrame.from_dict(dd)

    return res


def get_val(df, col, op, selvalue):
    """
    Function to estimate values

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.
    col : TYPE
        DESCRIPTION.
    op : TYPE
        DESCRIPTION.
    selvalue : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    idx = op(df[col], selvalue)
    sel = df[idx]

    return len(sel)/len(df)
