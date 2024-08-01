#!/usr/bin/env python
# coding: utf-8
#
## Plot WROMY data
##_______________________________

import os
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

if os.uname().nodename == 'lighthouse':
    root_path = '/home/andbro/'
    data_path = '/home/andbro/kilauea-data/'
    archive_path = '/home/andbro/freenas/'
    bay_path = '/home/andbro/bay200/'
elif os.uname().nodename == 'kilauea':
    root_path = '/home/brotzer/'
    data_path = '/import/kilauea-data/'
    archive_path = '/import/freenas-ffb-01-data/'
    bay_path = '/bay200/'
elif os.uname().nodename in ['lin-ffb-01', 'ambrym', 'hochfelln']:
    root_path = '/home/brotzer/'
    data_path = '/import/kilauea-data/'
    archive_path = '/import/freenas-ffb-01-data/'
    bay_path = '/bay200/'

##_______________________________
## Configurations

conf_path = root_path+'Documents/ROMY/HTML_Monitor/python_scripts/'

config_file = 'config.yaml'

# try:
#     config = __readYaml(conf_path, config_file)['WROMY']
# except:
#     print(f" -> failed to load configuration file!")

config = {}

config['channel'] = None

config['stations'] = [1, 4, 5, 6, 7, 8, 9]

config['colors'] = {  'WS1':'darkgreen', 'WS4':'purple', 'WS5':'darkred',
                        'WS6':'darkblue', 'WS7':'darkorange', 'WS8':'darkcyan', 'WS9':'cyan',}

config['tend'] = UTCDateTime.now()

config['timePeriod'] = 21 # days

config['tbeg'] = config['tend'] - config['timePeriod']*86400


config['path_to_wromy_data'] = archive_path+"romy_archive/"
config['path_to_furt_data'] = bay_path+'/gif_online/FURT/WETTER/'


# config['outpath'] = '/home/brotzer/'
# config['outname'] = 'html_wromy_plots'

# config['resample'] = 20



##_______________________________
## Methods

def __reply(msg):
    print(f"   -> {msg}")


def __read_wromy_data(config):
    '''
    reads data from T1 to T2
    '''
    path = f"{config['path_to_archive']}{config['tbeg'].year}/BW/WROMY/{config['channel']}.D/"

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



def __indicate_gaps_with_nan(df, config):

    differences = np.diff(df.totalSeconds, n=1)

    sample_time_errors = [j for j in differences if j != config['resample']]

    if len(sample_time_errors) != 0:
        print(f"  -> ERROR: Found {len(sample_time_errors)} errors for the sampling time!\n")

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


def __processing(data, config):

    filter_length = 10*config['resample']

    data.iloc[:,3:6] = data.iloc[:,3:6].rolling(filter_length).mean()
    __reply(f"Filter: rooling mean {filter_length}!")

    data = data[data.index % config['resample'] == 0]
    __reply(f"Resampling: keep every {config['resample']}nth sample!")

    data.fillna(np.nan)

    return data


def __read_furt_data(config):
    '''
    Load data of FURT wather station
    '''

    if not Path(config['path_to furt_data']).exists():
        print(f"  -> Path: {config['path_to_furt_data']}, does not exists!")
        return


    df = pd.DataFrame()

    for i, date in enumerate(np.arange(config['tbeg'].date, (config['tend']+86400).date)):

        date = UTCDateTime(str(date)).date
        filename = f'FURT.WSX.D.{str(date.day).rjust(2,"0")}{str(date.month).rjust(2,"0")}{str(date.year).rjust(2,"0")[-2:]}.0000'

        print(f'   reading {filename} ...')

        try:

            df0 = pd.read_csv(config['path_to_furt_data']+filename, usecols=[0,1,10,12,13], names=['date', 'time', 'T', 'H', 'P'])

            ## substitute strings with floats
            df0['T'] = df0['T'].str.split("=", expand=True)[1].str.split("C", expand=True)[0].astype(float)
            df0['P'] = df0['P'].str.split("=", expand=True)[1].str.split("H", expand=True)[0].astype(float)
            df0['H'] = df0['H'].str.split("=", expand=True)[1].str.split("P", expand=True)[0].astype(float)


            ## replace error ipath_to_archivendicating values (-9999, 999.9) with NaN values
            df0.replace(to_replace=-9999, value=nan, inplace=True)
            df0.replace(to_replace=999.9, value=nan, inplace=True)


            if df.empty:
                df = df0
            else:
                df = pd.concat([df, df0])
        except:
            print(f"  -> File: {filename}, does not exists!")

    df.reset_index(inplace=True, drop=True)

    return df



def __make_plot_all_stations_and_furt(data, furt, config, events=None):

    N = 3
    font = 13
    datasize = 0

    fig, axes = plt.subplots(N,1, figsize=[15,15], sharex=True)

    plt.subplots_adjust(hspace=0.1)

    max_val, min_val = np.zeros(N)*np.nan, np.zeros(N)*np.nan

    timeaxis_furt, ticks_furt, ticklabels_furt, text_furt = __get_timeaxis(
                                                                           dates=furt.iloc[:,0],
                                                                           times=furt.iloc[:,1],
                                                                           unit="date",
                                                                           unitmode="relative",
                                                                           dateformat="ddmmyy",
                                                                          )
    for station in data.keys():

        df = data.get(station)

        for u in range(3):
            maximum = df.iloc[:,u+3].dropna().max()
            minimum = df.iloc[:,u+3].dropna().min()
            if maximum > max_val[u] or np.isnan(max_val[u]):
                max_val[u] = maximum
            if minimum < min_val[u] or np.isnan(min_val[u]):
                min_val[u] = minimum

        timeaxis, ticks, ticklabels, text = __get_timeaxis(dates=df.iloc[:,1],
                                                           times=df.iloc[:,2],
                                                           unit="date",
                                                           unitmode="relative",
                                                           dateformat="yyyymmdd",
                                                          )



        ## select ticks and ticklabels for longest data series
        if df.shape[0] > datasize:
            datasize = df.shape[0]
            xticks = ticks
            xlabels = ticklabels
            timeaxis_min, timeaxis_max = timeaxis[0], timeaxis[-1]

        ## plot data and adjust axes automatically
        for i in range(N):
#             axes[i].scatter(timeaxis, df.iloc[:,i+3], s=1, color='grey',lw=3, zorder=2)
            axes[i].plot(timeaxis, df.iloc[:,i+3], color=config['colors'][station], lw=1.5, zorder=2, label=station)

            if station == list(data.keys())[-1]:
                axes[i].plot(timeaxis_furt, furt.iloc[:,i+2], color='darkgrey', lw=1.5, zorder=1, label="FURT")

            axes[i].grid(ls="--",color='grey', zorder=0)

            if i == 0:
                axes[i].set_ylabel("Temperature (°C)",fontsize=font)
            elif i == 1:
                axes[i].set_ylabel("Air Pressure (hPa)",fontsize=font)
            elif i == 2:
                axes[i].set_ylabel("rel. Humidity (%)",fontsize=font)
                axes[i].set_xlim(timeaxis_min, timeaxis_max)

            if events:
                for event in events:
#                     axes[i].axvline(event, color='r', zorder=0, ls="-.")
                    axes[i].axvspan(event[0], event[1], color="lightgrey", alpha=0.4, zorder=1)

        axes[N-1].set_xticklabels(xlabels)
        axes[N-1].set_xticks(xticks)
        axes[N-1].set_xlabel(text, fontsize=font)
        axes[N-1].legend(loc='upper center', ncol=7+1, bbox_to_anchor=(0.5, -0.15), fancybox=True)

    #plt.show();
    return fig

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Main
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if __name__ == "__main__":

    ## _______________________________
    ## Loading

    furt = __read_furt_data(config)
    print(furt)
    furt = furt[['date', 'time', 'T', 'P', 'H']]
    furt = __processing(furt, config)

    data = {}

    for i in config['stations']:

        config['channel'] = 'WS'+str(i)

        ## load data as DataFrame
        df_new = __read_wromy_data(config)

        ## check for gaps
        # df_new = __indicate_gaps_with_nan(df_new, config)

        ## processing
        df_new = __processing(df_new, config)

        print(df_new.head())
        ## add to dictionary
        data[config.get('channel')] = df_new; del df_new


    ## _______________________________
    ## Plotting

    fig = __make_plot_all_stations_and_furt(data, furt, config);

    ## _______________________________
    ## Saving

    __savefig(fig, outpath=config['outpath'], outname=config['outname'], mode='png', dpi=500)

## END OF FILE