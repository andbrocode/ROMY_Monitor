#!/usr/bin/env python
# coding: utf-8

# ## Plot all traces of ROMY array


##_____________________________________________________________
## Import

import obspy as obs
import numpy as np
import matplotlib.pyplot as plt

from time import sleep
from numpy import isnan, ones, nan

from andbro__querrySeismoData import __querrySeismoData
from andbro__get_timeaxis import __get_timeaxis
from andbro__savefig import __savefig

import warnings
warnings.filterwarnings("ignore")

##_____________________________________________________________
## Methods

def __get_stream(config):

    st = obs.Stream()

    print("")

    for i in range(len(config.get("sta"))):

        seed = f'{config.get("net")[i]}.{config.get("sta")[i]}.{config.get("loc")}.{config.get("cha")}'

        print(f"loading {seed}...")

        try:
            st0, inv = __querrySeismoData(
                                        seed_id=seed,
                                        starttime=config.get("tbeg"),
                                        endtime=config.get("tend"),
                                        where=config.get("repository"),
                                        path=None,
                                        restitute=True,
                                        detail=None,
                                        fill_value=None,
                                        )

            if len(st0) != 0:
                st += st0
        except Exception as e:
            print("  -> Failed!")
            print(e)

    return st


def __empty_trace(config, i):
    tr = obs.Trace()
    tr.stats.sampling_rate = 20
    tr.data = np.zeros(int((config['tend']-config['tbeg'])*tr.stats.sampling_rate)+1) * np.nan
    tr.stats.station = config['sta'][i]
    tr.stats.network = config['net'][i]
    tr.stats.channel = config['cha']
    tr.stats.starttime = config['tbeg']

    return tr

def __makeplot_array_traces(st, config):
    switch = True

    N = len(config['sta'])

    fig, axes = plt.subplots(N, 1, figsize=(15,20), sharex='col')

    plt.subplots_adjust(hspace=0.2,wspace=0.2)

    for i, sta in enumerate(config['sta']):

        try:
            tr = st.select(station=sta)[0]
        except:
            tr = __empty_trace(config, i)



        print(f' plotting {tr.stats.station} ...')
        try:
            timeaxis, ticks, ticklabels, text = __get_timeaxis(
                                                             utcdatetime=tr.times(type="utcdatetime"),
                                                             unit="time",
                                                             unitmode="absolute",
                                                             dateformat="yyyymmdd",
                                                              )

            sta_cha = f"{tr.stats.station}.{tr.stats.channel}"

            axes[i].plot(timeaxis, tr.data, color='k', label=sta_cha, lw=0.6, zorder=2)

            if switch:
                switch = not switch
                tmin, tmax = timeaxis[0], timeaxis[-1]
            if tmax < timeaxis[-1]:
                tmax = timeaxis[-1]
        except:
            print("nothing to plot")



        axes[i].legend(loc="upper right")

        if i == N-1:
            axes[i].set_xticks(ticks)
            axes[i].set_xticklabels(ticklabels)
            axes[i].set_xlabel(text)
            axes[i].set_xlim(tmin, tmax)

    plt.show()
    return fig

def __processing(st, config):

    ## check for NaN in data and create masks
    masks, masks_empty = [], True
    for tr in st:
        if isnan(tr.data).any():
            mask = ones(len(tr.data))
            for i, e in enumerate(tr.data):
                if isnan(e):
                    tr.data[i] = 0
                    mask[i] = nan
            print("created masks for NaN values ")
            masks.append(mask)
            masks_empty = False
        else:
            masks.append([])

        ## Filtering
        if config['setFilter']:
            if config.get("filter_type") in ['bp', 'bandpass']:
                tr=tr.filter("bandpass", freqmin=config.get("filter_corners")[0], freqmax=config.get("filter_corners")[1], corners=4, zerophase=True)
            elif config.get("filter_type") in ['lp', 'lowpass']:
                tr=tr.filter("lowpass", freq=config.get("filter_corners")[1], corners=4, zerophase=True)
            elif config.get("filter_type") in ['hp', 'highpass']:
                tr=tr.filter("highpass", freq=config.get("filter_corners")[0], corners=4, zerophase=True)
        ## reapply masks
        if not masks_empty:
            for i, tr in enumerate(st):
                if i < len(masks):
                    if not len(masks[i]) == 0:
                        tr.data *= masks[i]

    return st

def __update_traces(st, config):

    config['tbeg_old'], config['tend_old'] = config['tbeg'], config['tend']
    config['tbeg'] = config['tend']
    config['tend'] = obs.UTCDateTime.now() - config['lagTime']

    st_update = __get_stream(config)

    if len(st_update) > 0:
        st += st_update

    st = st.merge(method=1, fill_value=nan)

    config['tbeg'] = config['tbeg_old'] + (config['tend'] - config['tbeg'])

    st = st.trim(config['tbeg'], config['tend'], pad=True, fill_value=nan)

    return st, config


##_____________________________________________________________
## Configurations

config = {}

config['lagTime'] = 30*60  ## seconds
config['showTime'] = 30*60 ## seconds

config['tbeg'] = obs.UTCDateTime().now() - config['lagTime'] - config['showTime']
config['tend'] = obs.UTCDateTime().now() - config['lagTime']

config['sta'] = ['GELB','GRMB','BIB','TON', 'ALFT', 'FFB1', 'FFB2', 'FFB3', 'FUR']
config['net'] = ['BW','BW','BW','BW','BW','BW','BW','BW','GR']
config['loc'] = ''
config['cha'] = 'BHZ'

config['repository'] = "jane"


config['setFilter'] = True
config['filter_type'] = 'bandpass'
config['filter_corners'] = [0.01, 0.3]

config['outpath'] = '/home/brotzer/Documents/ROMY/HTML_Monitor/figures/'
config['outname'] = 'test'

##_____________________________________________________________
## Load Data

st = __get_stream(config)

i = 0
while True:

    if i > 0:
        st, config = __update_traces(st, config)
    i+=1

    st = __processing(st, config)

    fig = __makeplot_array_traces(st, config);

    print("\nwaiting ... \n"); sleep(120)

    ## ______________________________________
    ## Saving

    __savefig(fig, outpath=config['outpath'], outname=config['outname'], mode='png', dpi=300)

##_____________________________________________________________
## END OF FILE