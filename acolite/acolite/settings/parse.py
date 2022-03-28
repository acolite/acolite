## def parse
## parses acolite settings files and merges with defaults
## written by Quinten Vanhellemont, RBINS
## 2021-03-09
## modifications: 2022-03-28 (QV) moved int and float lists to external files


def parse(sensor, settings=None, merge=True):
    import os, time
    import numpy as np
    import scipy.ndimage
    import acolite as ac

    ## read default settings for sensor
    if (sensor is not None) | (merge):
        setu = ac.acolite.settings.load(sensor)
    else:
        setu = {}

    ## add user settings
    if settings is not None:
        if type(settings) is str:
            sets = ac.acolite.settings.read(settings)
            if merge:
                for k in sets: setu[k] = sets[k]
            else:
                setu = {k:sets[k] for k in sets}
        elif type(settings) is dict:
            if merge:
                for k in settings: setu[k] = settings[k]
            else:
                setu = {k:settings[k] for k in settings}

    ## make sure luts setting is a list
    if 'luts' in setu:
        if type(setu['luts']) is not list: setu['luts'] = [setu['luts']]

    ## import settings that need to be converted to ints and floats
    int_list = ac.acolite.settings.read_list(ac.config['data_dir']+'/ACOLITE/settings_int.txt')
    float_list = ac.acolite.settings.read_list(ac.config['data_dir']+'/ACOLITE/settings_float.txt')

    ## convert values to numbers
    for k in setu:
        if k not in setu: continue
        if setu[k] is None: continue

        if type(setu[k]) is list:
            if k in int_list: setu[k] = [int(i) for i in setu[k]]
            if k in float_list: setu[k] = [float(i) for i in setu[k]]
        else:
            if k in int_list: setu[k] = int(setu[k])
            if k in float_list: setu[k] = float(setu[k])

    ## default pressure
    if 'pressure' in setu:
        setu['pressure'] = 1013.25 if setu['pressure'] is None else float(setu['pressure'])

    return(setu)
