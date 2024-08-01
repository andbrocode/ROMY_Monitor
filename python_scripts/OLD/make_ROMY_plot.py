# -*- coding: utf-8 -*-
"""
Spyder Editor

Plots all ROMY Z,U,V,W traces + RLAS Z

"""
__author__ = 'AndreasBrotzer'

##_____________________________________________________________

'''---- import libraries ----'''

import obspy as obs
import sys
import matplotlib.pyplot as plt

from numpy import floor, isnan, nanmean, nanmax, nanmin, ones, nan
from pathlib import Path
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

config = __readYaml(conf_path, conf_name)['ROMY']


config['tend'] = UTCDateTime.now()
config['tbeg'] = config['tend'] - config['timePeriod']*60


##_____________________________________________________________
'''---- define methods ----'''


def __requesting_data(user_params):

    print("\n"+reply+f"requesting data for {user_params['seed_id']} from {user_params['repository']}...\n")

    st = obs.Stream()
    inv = []
    try:

        for seed_id in user_params.get('seed_id'):
            st0, inv0 = __querrySeismoData(
                                        seed_id=seed_id,
                                        starttime=user_params.get("tbeg"),
                                        endtime=user_params.get("tend"),
                                        where=user_params.get("repository"),
                                        path=user_params.get("datapath"),
                                        restitute=True,
                                        detail=None,
                                        )
            st += st0
            inv.append(inv0)

        if len(st) != 0:
            print(reply + "loaded data")
        else:
            print(reply + "empty stream!")
            sys.exit()

    except Exception as e:
        print(e)
        sys.exit()

    else:
        for tr in st:
            if isnan(tr.data).any():
                tr.data -= nanmean(tr.data)
            else:
                tr.detrend("demean")
    # main()
        print(reply + "applied demean ")

    return st, inv


def __makeplot_romy_traces(st, config):


    N = 5
    font = 13

    fig, axes = plt.subplots(N, 1, figsize=(15,10), sharex='col')

    plt.subplots_adjust(hspace=0.2,wspace=0.2)

    ## check for appropriate time axis unit
    time_factor = {'sec':1, 'min':60, 'hr':3600}
    time_unit = 'sec'
    if st[0].stats.npts * st[0].stats.delta > 3600:
        time_unit = 'min'
    elif st[0].stats.npts * st[0].stats.delta > 12*3600:
        time_unit = 'hr'

    ## _______________________________________________

    for i, tr in enumerate(st):

        timeaxis, ticks, ticklabels, text = __get_timeaxis(
                                                         utcdatetime=tr.times(type="utcdatetime"),
                                                         unit="time",
                                                         unitmode="absolute",
                                                         dateformat="yyyymmdd",
                                                        )

        sta_cha = f"{tr.stats.station}.{tr.stats.channel}"

        # axes[i].plot(tr.times()/time_factor.get(time_unit),tr.data, color='k', label=sta_cha, lw=0.6, zorder=2)
        axes[i].plot(timeaxis,tr.data, color='k', label=sta_cha, lw=0.9, zorder=2)

        gaps = isnan(tr)
        bounds = [0, 1]
        if not isnan(nanmin(tr.data)):
            bounds[0] = nanmin(tr.data)
        if not isnan(nanmax(tr.data)):
            bounds[1] = nanmax(tr.data)

        # axes[i].fill_between(tr.times()/time_factor.get(time_unit), bounds[0], bounds[1], where=gaps, alpha=0.4, color='red', zorder=1)
        axes[i].fill_between(timeaxis, bounds[0], bounds[1], where=gaps, alpha=0.4, color='red', zorder=1)

        axes[i].set_ylabel(r'$\Omega$ (rad/s)', fontsize=font-2)

        axes[i].legend(loc="upper right", fontsize=font, bbox_to_anchor=(1.05, 1.20), shadow=True)
        # axes[i].grid(color="grey", zorder=0, ls="--")

        if config.get('freq'):
            if len(config['freq']) == 1:
                axes[0].set_title(f"lowpass: {config['freq'][0]} Hz ")
            elif len(config['freq']) == 2:
                axes[0].set_title(f"bandpass: {config['freq'][0]} - {config['freq'][1]} Hz ")


        if i == N-1:
            # axes[i].set_xlabel(f"Time from {str(config['tbeg'])[0:10]} {str(config['tbeg'])[11:16]} UTC ({time_unit})")
            axes[i].set_xticks(ticks)
            axes[i].set_xticklabels(ticklabels, fontsize=font-1)
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
                tr = tr.filter("bandpass", freqmin=config.get("filter_corners")[0], freqmax=config.get("filter_corners")[1], corners=4, zerophase=True)
            elif config.get("filter_type") in ['lp', 'lowpass']:
                tr = tr.filter("lowpass", freq=config.get("filter_corners")[1], corners=4, zerophase=True)
            elif config.get("filter_type") in ['hp', 'highpass']:
                tr = tr.filter("highpass", freq=config.get("filter_corners")[0], corners=4, zerophase=True)

        ## reapply masks
        if not masks_empty:
            for i, tr in enumerate(st):
                if not len(masks[i]) == 0:
                    tr.data *= masks[i]

    return st

##_____________________________________________________________

if __name__ == "__main__":

    # config = __user_interaction()

    st, inv = __requesting_data(config)

    ## ______________________________________
    ## Processing

    st = __processing(st, config)

    ## ______________________________________
    ## Plotting

    fig = __makeplot_romy_traces(st, config)

    ## ______________________________________
    ## Saving

    __savefig(fig, outpath=config['outpath'], outname=config['outname'], mode='png', dpi=300)


##_____________________________________________________________
# End of File

