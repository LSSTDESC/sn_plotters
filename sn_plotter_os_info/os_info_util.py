from . import plt
import pandas as pd
import numpy as np
import glob

def plot_summary(data, field='DD:COSMOS',
                 varx='year', labx='year',
                 vary='seq_tot_y', laby='',
                 figtit='DD:COSMOS',
                 df_config=pd.DataFrame(),
                 cumsum=False):
    """
    Summary plot

    Parameters
    ----------
    data : pandas df
        Data to process.
    field : str, optional
        Field type. The default is 'DD:COSMOS'.
    varx : str, optional
        x-axis variable. The default is 'year'.
    labx : str, optional
        x-axis label. The default is 'year'.
    vary : str, optional
        y-axis variable. The default is 'seq_tot_y'.
    laby : str, optional
        y-axis label. The default is ''.
    figtit : str, optional
        Figure title. The default is 'DD:COSMOS'.
    df_config : pandas df, optional
        config for the plot. The default is pd.DataFrame().
    cumsum: bool, optional.
        to plot cumsum. The default is False.

    Returns
    -------
    None.

    """

    fig, ax = plt.subplots(figsize=(16, 8))
    fig.suptitle(figtit)
    fig.subplots_adjust(right=0.75)

    idx = data['target_name'] == field

    sel = data[idx]

    dbNames = sel['dbName'].unique()

    for dbName in dbNames:
        io = sel['dbName'] == dbName
        selb = sel[io]
        idxb = df_config['dbName'] == dbName
        selp = df_config[idxb]
        ls = selp['ls'].values[0]
        marker = selp['marker'].values[0]
        color = selp['color'].values[0]
        dbNameb = selp['dbName_plot'].values[0]

        tp = selb[vary]
        if cumsum:
            tp = np.cumsum(selb[vary])

        ax.plot(selb[varx], tp,
                ls=ls, marker=marker, color=color, mfc='None', label=dbNameb)

    ax.grid(visible=True)
    ax.set_xlabel(r'{}'.format(labx))
    ax.set_ylabel(r'{}'.format(laby))
    if laby == '':
        ax.tick_params(axis='y', labelrotation=20, labelsize=10)

    ax.legend(loc='upper center',
              bbox_to_anchor=(1.20, 0.7),
              ncol=1, fontsize=12, frameon=False)

    # plt.tight_layout()
    return ax

def load_data(dbDir, dbName, fields, fieldType='DD'):
    """
    Function to load the data

    Parameters
    ----------
    dbDir : str
        data main directory.
    dbName : str
        OS to analyze.
    fields : list(str)
        list of fields to process.
    fieldType : str, optional
        type of field (DD/WFD). The default is 'DD'.

    Returns
    -------
    df : pandas df
        output data.

    """

    theDir = '{}/{}'.format(dbDir, dbName)
    prefix = '{}_pixels_{}'.format(fieldType, dbName)
    
    df = pd.DataFrame()
    fName = '{}/{}*.hdf5'.format(theDir, prefix)
    print('search_path',fName)
    fis = glob.glob(fName)
    for fi in fis:
        dd = pd.read_hdf(fi)
        df = pd.concat((df, dd))
    
    return df

def load_sel_dust(dbDir,dbNames,fields,fieldType='DD',
              nside=128,ebvofMW_max=0.25):
    
    df = pd.DataFrame()
    for dbName in dbNames:
        df_ = load_data(dbDir, dbName, fields, fieldType=fieldType)
        df_['dbName'] = dbName
        df_['nvisits_gri'] = df_['nvisits_g']+df_['nvisits_r']+df_['nvisits_i']
        df = pd.concat((df,df_))

    #load the dust_map

    fDust = 'reference_files/dustmap_{}_delta_mag_dust.hdf5'.format(nside)
    df_dust = pd.read_hdf(fDust)

    df = df.merge(df_dust,left_on=['healpixID'],right_on=['healpixID'])

    #select on dust

    idx = df['ebvofMW'] < ebvofMW_max
    df = df[idx]

    print('nentries',len(df))
    
    return df