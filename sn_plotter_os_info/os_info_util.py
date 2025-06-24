from . import plt
import pandas as pd


def plot_summary(data, field='DD:COSMOS',
                 varx='year', labx='year',
                 vary='seq_tot_y', laby='',
                 figtit='DD:COSMOS',
                 df_config=pd.DataFrame()):
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

        ax.plot(selb[varx], selb[vary],
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
