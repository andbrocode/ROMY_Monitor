# -*- coding: utf-8 -*-

"""
Spyder Editor

This script enables to plot a period of the tiltmeter data located at ROMY.

"""
__author__ = 'AndreasBrotzer'
__year__   = '2021'


##_____________________________________________________________
'''---- import libraries ----'''

import matplotlib.pyplot as plt
import obspy as obs
import sys
import yaml

from numpy import floor, isnan, nanmean, nanmax, nanmin, ones, nan, array
from obspy import UTCDateTime

from andbro__querrySeismoData import __querrySeismoData
from andbro__savefig import __savefig
from andbro__readYaml import __readYaml

import warnings
warnings.filterwarnings("ignore")

reply = "\n    -> "

##____________________________________________________________
'''---- Configurations ----'''

## configurations for the plot
conf_path = '/home/brotzer/Documents/ROMY/HTML_Monitor/python_scripts/'
conf_name = 'config.yaml'

config = __readYaml(conf_path, conf_name)['TILT']


config['tend'] = UTCDateTime.now()
config['tbeg'] = config['tend'] - config['timePeriod']*86400

## specifications of Tiltmeters
confTILT_path = ''
confTILT_name = ''

##_____________________________________________________________
'''---- Methods ----'''


###########################################################
## Loading

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
                                        fill_value=-9999,
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

    # else:
    #     for tr in st:
    #         if isnan(tr.data).any():
    #             tr.data -= nanmean(tr.data)
    #         else:
    #             tr.detrend("demean")
    #     print(reply + "applied demean ")

    return st, inv

###########################################################
## Plotting

def __makeplot_tilt_traces(st, config):

    N = 6
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

        sta_cha = f"{tr.stats.station}.{tr.stats.channel}"
        if tr.stats.channel[0] == "M":
            name = 'Tilt PT'
        elif tr.stats.channel[0] == "L":
            name = 'Tilt BT'
        if tr.stats.channel[-1] in ["N","E"]:
            comp, unit = tr.stats.channel[-1], '(rad)'
        else:
            comp, unit = 'T', '(°C)'

        axes[i].plot(tr.times()/time_factor.get(time_unit),tr.data, color='k', label=sta_cha, lw=0.6, zorder=2)

        gaps = isnan(tr)
        bounds = [0, 1]
        if not isnan(nanmin(tr.data)):
            bounds[0] = nanmin(tr.data)
        if not isnan(nanmax(tr.data)):
            bounds[1] = nanmax(tr.data)

        axes[i].fill_between(tr.times()/time_factor.get(time_unit), bounds[0], bounds[1], where=gaps, alpha=0.4, color='red', zorder=1)

        axes[i].set_ylabel(f'{name} {comp} {unit}')

        axes[i].legend(loc="upper right")
        axes[i].grid(color="grey", zorder=0, ls="--")

        if config.get('freq'):
            if len(config['freq']) == 1:
                axes[0].set_title(f"lowpass: {config['freq'][0]} Hz ")
            elif len(config['freq']) == 2:
                axes[0].set_title(f"bandpass: {config['freq'][0]} - {config['freq'][1]} Hz ")


        if i == N-1:
            axes[i].set_xlabel(f"Time from {str(config['tbeg'])[0:10]} {str(config['tbeg'])[11:16]} UTC ({time_unit})")
            axes[i].set_xlim(tr.times()[0]/time_factor.get(time_unit), tr.times()[-1]/time_factor.get(time_unit))

    #plt.show()

    return fig

###########################################################
## Processing

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


        ## ______________________________________
        ## Filtering

    for tr in st:
        if tr.stats.channel[-1] != 'T':
            tr.detrend('demean')

            if config['setFilter']:
                tr.taper(0.1)

                if config.get("filter_type") in ['bp', 'bandpass']:
                    tr.filter("bandpass", freqmin=config.get("filter_corners")[0], freqmax=config.get("filter_corners")[1], corners=4, zerophase=True)
                elif config.get("filter_type") in ['lp', 'lowpass']:
                    tr.filter("lowpass", freq=config.get("filter_corners")[1], corners=4, zerophase=True)
                elif config.get("filter_type") in ['hp', 'highpass']:
                    tr.filter("highpass", freq=config.get("filter_corners")[0], corners=4, zerophase=True)
        ## ______________________________________
        ## Resampling

    for tr in st:

        if config.get("resample") in ['yes','y']:
            if tr.stats.channel[-1] != 'T':
    #            tr.detrend('demean')
                tr.resample(config.get("resampling_frequency"), window='hanning')
            print(reply+ f"resampled data with {config.get('resampling_frequency')} Hz")


    ## reapply masks
    if not masks_empty:
        for i, tr in enumerate(st):
            if not len(masks[i]) == 0:
                tr.data *= masks[i]

    print(reply+ "processed data (demean + filter)")
    return st

###########################################################
## Conversion

def __conversion(st, confBT, confPT):

    def convertTemp(trace, gain):
        Tvolt = trace.data * gain
        return 10.18 - 11.59* Tvolt + 0.3335* Tvolt**2 - 0.5316* Tvolt**3

    def convertTilt(trace, conversion, sensitivity):
        return trace.data * conversion * sensitivity
        # print( type(conversion), type(sensitivity), type(trace.data) )

    for tr in st:
        if tr.stats.channel == 'MAT':
            tr.data = convertTemp(tr, confPT['gainTemp'])
        elif tr.stats.channel == 'MAN':
            tr.data = convertTilt(tr, confPT['convPTN'], confPT['gainTilt'])
        elif tr.stats.channel == 'MAE':
            tr.data = convertTilt(tr, confPT['convPTE'], confPT['gainTilt'])

        elif tr.stats.channel == 'LAT':
            tr.data = convertTemp(tr, confBT['gainTemp'])
        elif tr.stats.channel == 'LAN':
            tr.data = convertTilt(tr, confBT['convBTN'], confBT['gainTilt'])
        elif tr.stats.channel == 'LAE':
            tr.data = convertTilt(tr, confBT['convBTE'], confBT['gainTilt'])

    print(reply + "converted data")
    return st

###########################################################
## Main


##_____________________________________________________________
if __name__ == "__main__":

    # config = __user_interaction()

    ## specification of the tiltmeters
    confTilt = __readYaml('/home/brotzer/Documents/ROMY/','tiltmeter.conf')
    confPT = confTilt['PT']
    confBT = confTilt['BT']

    ## ______________________________________
    ## Request Data

    st, inv = __requesting_data(config)

    ## ______________________________________
    ## Conversion

    st = __conversion(st, confBT, confPT)

    ## ______________________________________
    ## Processing
    if not config['setFilter'] in ['n', 'no', None]:
        st = __processing(st, config)

    ## ______________________________________
    ## Plotting

    fig = __makeplot_tilt_traces(st, config)


    ## ______________________________________
    ## Saving

    __savefig(fig, outpath=config['outpath'], outname=config['outname'], mode='png', dpi=500)



##_____________________________________________________________
# End of File