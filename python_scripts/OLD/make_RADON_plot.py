#!/usr/bin/env python
# coding: utf-8
#
# Plot RADON data of WS1
#
## __________________________________________
## Importing

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from numpy import nan
from pathlib import Path
from obspy import UTCDateTime

from andbro__get_timeaxis import __get_timeaxis
from andbro__savefig import __savefig
from andbro__readYaml import __readYaml

import warnings
warnings.filterwarnings("ignore")

## __________________________________________
## Configurations

conf_path = '/home/brotzer/Documents/ROMY/HTML_Monitor/python_scripts/'
conf_name = 'config.yaml'

config = __readYaml(conf_path, conf_name)['RADON']



# config = {}

# config['channel'] = 'RDN'

config['tend'] = UTCDateTime.now()
config['tbeg'] = config['tend'] - config['timePeriod']*86400

# config['pathToData'] = f'/import/freenas-ffb-01-data/romy_archive/'

# config['resample'] = 4


## __________________________________________
## Methods

## ##########################################################
def __reply(msg):
    print(f"   -> {msg}")

## ##########################################################
def __read_wromy_data(config):
    '''
    reads data from T1 to T2
    '''

    path = f"{config['pathToData']}{config['tbeg'].year}/BW/WROMY/{config['channel']}.D/"

    if not Path(path).exists():
        __reply(f"Path: {path}, does not exists!")
        return

    j1, j2  = config['tbeg'].julday, config['tend'].julday
    year    = config['tbeg'].year
    cha     = config['channel']
    df = pd.DataFrame()

    for doy in range(j1, j2+1):

        fileName = f'BW.WROMY.{cha}.D.{year}.{doy}'

        print(f'   reading {fileName} ...')

        try:
            df0 = pd.read_csv(path+fileName)

            ## replace error indicating values (-9999, 999.9) with NaN values
            df0.replace(to_replace=-9999, value=nan, inplace=True)
            df0.replace(to_replace=999.9, value=nan, inplace=True)


            if doy == j1:
                df = df0
            else:
                df = pd.concat([df,df0])
        except:
            __reply(f"File: {fileName}, does not exists!")

    df.reset_index(inplace=True, drop=True)

    ## add columns with total seconds
    if 'Seconds' in df.columns:
        totalSeconds = df.Seconds + (df.Date - df.Date.iloc[0]) * 86400
        df['totalSeconds'] = totalSeconds


    __reply("Done \n")

    return df

## ##########################################################
def __indicate_gaps_with_nan(df, config):

    differences = np.diff(df.totalSeconds, n=1)


    ## ______________

    sample_time_errors = [j for j in differences if j != config['resample']]

    if len(sample_time_errors) != 0:
        print(f"  -> ERROR: Found {len(sample_time_errors)} errors for the sampling time!\n")


    ## ______________

    gaps = [list(differences).index(k) for k in differences if k > 2*config['resample']] or []
    if gaps and gaps[0] in [0, 0.0]:
        gaps.pop(0)
    del differences

    for x in gaps:
        fill_row = [i+config['resample'] if n not in [3,4,5] else np.nan for n, i in enumerate(df.iloc[x,:])]
        fill_row[0] = int(df.iloc[x,0])
        fill_row[1] = int(df.iloc[x,1])
        fill_row[2] = int(df.iloc[x,2])
        df.loc[x+0.5] = fill_row


    df = df.sort_index().reset_index(drop=True).convert_dtypes()

    print(f"  -> Marked {len(gaps)} gaps with NaN values!\n")

    return df


## ##########################################################
def __processing(data, config):

    filter_length = 10*config['resample']

    data.iloc[:,3:6] = data.iloc[:,3:6].rolling(filter_length).mean()
    __reply(f"Filter: rooling mean {filter_length}!")

    data = data[data.index % config['resample'] == 0]
    __reply(f"Resampling: keep every {config['resample']}nth sample!")

    return data


## ##########################################################
def __makeplot_radon(data, events=None):

    data_to_plot = [5,4,3]

    N = len(data_to_plot)
    font = 13
    datasize = 0


    timeaxis, ticks, ticklabels, text = __get_timeaxis(
                                                       dates=data.iloc[:,1],
                                                       times=data.iloc[:,2],
                                                       unit="date",
                                                       unitmode="relative",
                                                       dateformat="yyyymmdd",
                                                      )


    fig, axes = plt.subplots(N,1, figsize=(15,8), sharex=True)

    for i in range(N):

        axes[i].plot(timeaxis, data.iloc[:,data_to_plot[i]], label=data.columns[data_to_plot[i]], color='k', zorder=2)

        axes[i].grid(ls="--", color='grey', zorder=0, alpha=0.7)

        axes[i].set_ylabel(data.columns[data_to_plot[i]], fontsize=font)

        axes[i].set_ylim(0, np.nanmax(data.iloc[:,data_to_plot[i]])+20)

        if i in range(2):
            axes[i].axhspan(0,100, color='green', alpha=0.3, zorder=0)
            axes[i].axhspan(100,150, color='yellow', alpha=0.3, zorder=0)
            axes[i].axhspan(150,np.nanmax(data.iloc[:,data_to_plot[i]])+20, color='red', alpha=0.3, zorder=0)

            axes[i].axhline(300, color='darkorange', alpha=0.9, zorder=1)

        if events:
            for event in events:
                if UTCDateTime(event[0]) > timeaxis[0]:
                    axes[i].axvspan(event[0], event[1], color="grey", alpha=0.4, zorder=0)

    axes[N-1].set_xticks(ticks)
    axes[N-1].set_xticklabels(ticklabels)
    axes[N-1].set_xlabel(text, fontsize=font)

    return fig


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Main
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if __name__ == "__main__":

    ## _______________________________
    ## loading

    radon = __read_wromy_data(config)

    # radon = __indicate_gaps_with_nan(radon, config)

    radon.insert(0, 'Seconds', np.zeros(radon.shape[0]))

    radon = __processing(radon, config)

    ## _______________________________
    ## Plotting

    fig = __makeplot_radon(radon, events=None);

    ## _______________________________
    ## Saving

    __savefig(fig, outpath=config['outpath'], outname=config['outname'], mode='png', dpi=500)



## END OF FILE