# -*- coding: utf-8 -*-
"""
Spyder Editor

Plots all traces of ROMY Array

"""
__author__ = 'AndreasBrotzer'

##_____________________________________________________________

'''---- import libraries ----'''

import obspy as obs
import matplotlib.pyplot as plt

from numpy import floor, isnan, nanmean, nanmax, nanmin, ones, nan, zeros
from obspy import UTCDateTime

from andbro__querrySeismoData import __querrySeismoData
from andbro__savefig import __savefig
from andbro__readYaml import __readYaml
from andbro__get_timeaxis import __get_timeaxis

import warnings
warnings.filterwarnings("ignore")

reply = "\n    -> "


##____________________________________________________________
'''---- Configurations ----'''

## configurations for the plot
conf_path = '/home/brotzer/Documents/ROMY/HTML_Monitor/python_scripts/'
conf_name = 'config.yaml'

config = __readYaml(conf_path, conf_name)['ARRAY']


config['tbeg'] = UTCDateTime().now() - 3600
config['tend'] = UTCDateTime().now()


##_____________________________________________________________
'''---- define methods ----'''

def __get_stream(config):
    st = obs.Stream()

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
        except:
            print("  -> Failed!")
    return st

def __empty_trace(config, i):
    tr = obs.Trace()
    tr.stats.sampling_rate = 20
    tr.data = zeros(int((config['tend']-config['tbeg'])*tr.stats.sampling_rate)+1) * nan
    tr.stats.station = config['sta'][i]
    tr.stats.network = config['net'][i]
    tr.stats.channel = config['cha']
    tr.stats.starttime = config['tbeg']

    return tr

def __makeplot_array_traces(st, config):

    N = len(config['sta'])

    fig, axes = plt.subplots(N, 1, figsize=(18,20), sharex='col')

    plt.subplots_adjust(hspace=0.2,wspace=0.2)

    for i, sta in enumerate(config['sta']):

        try:
            tr = st.select(station=sta)[0]
        except:
            tr = __empty_trace(config, i)

        print(f'plotting {tr.stats.station} ...')

        timeaxis, ticks, ticklabels, text = __get_timeaxis(
                                                         utcdatetime=tr.times(type="utcdatetime"),
                                                         unit="time",
                                                         unitmode="absolute",
                                                         dateformat="yyyymmdd",
                                                          )

        sta_cha = f"{tr.stats.station}.{tr.stats.channel}"

        axes[i].plot(timeaxis, tr.data, color='k', label=sta_cha, lw=0.9, zorder=2)

        axes[i].legend(loc="upper right")

        if i == N-1:
            axes[i].set_xticks(ticks)
            axes[i].set_xticklabels(ticklabels)
            axes[i].set_xlabel(text)
            axes[i].set_xlim(timeaxis[0], timeaxis[-1])

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
            print(reply+"created masks for NaN values ")
            masks.append(mask)
            masks_empty = False
        else:
            masks.append([])

        ## Filtering
        if config['setFilter']:
            if config.get("filter_type") in ['bp', 'bandpass']:
                tr.filter("bandpass", freqmin=config.get("filter_corners")[0], freqmax=config.get("filter_corners")[1], corners=4, zerophase=True)
            elif config.get("filter_type") in ['lp', 'lowpass']:
                tr.filter("lowpass", freq=config.get("filter_corners")[1], corners=4, zerophase=True)
            elif config.get("filter_type") in ['hp', 'highpass']:
                tr.filter("highpass", freq=config.get("filter_corners")[0], corners=4, zerophase=True)

        ## reapply masks
        if not masks_empty:
            for i, tr in enumerate(st):
                    if len(masks) > 0:
                        if not len(masks[i]) == 0:
                            tr.data *= masks[i]


    return st

##_____________________________________________________________

if __name__ == "__main__":

    ## ______________________________________
    ## Loading

    st = __get_stream(config)

    ## ______________________________________
    ## Processing

    st = __processing(st, config)

    ## ______________________________________
    ## Plotting

    fig = __makeplot_array_traces(st, config);

    ## ______________________________________
    ## Saving

    __savefig(fig, outpath=config['outpath'], outname=config['outname'], mode='png', dpi=300)


##_____________________________________________________________
# End of File