## def acolite_settings
##
##
## Written by Quinten Vanhellemont 2017-11-30
## Last modifications:
##                2018-07-18 (QV) changed acolite import name

def load(settings):
    import os, glob
    import acolite as ac

    ## read defaults
    default_settings = '{}/config/defaults.txt'.format(ac.path)
    setd = ac.acolite.settings.read(default_settings)

    ## read settings file
    if settings is not None:
        ## path to settings file given
        if type(settings) is str:
            setf = '{}/config/defaults/{}.txt'.format(ac.path,settings)
            if (os.path.exists(settings)) and (not os.path.isdir(settings)):
                setu = ac.acolite.settings.read(settings)
            elif os.path.exists(setf):
                setu = ac.acolite.settings.read(setf)
            else:
                print('Settings file {} not found.'.format(settings))
                setu = setd
        elif type(settings) is dict:
            setu = settings
        else:
            print('Settings not recognised.')
            setu = setd
    else: setu={}

    ## set defaults for options not specified
    for key in setd.keys():
        if key in setu.keys(): continue
        else: setu[key] = setd[key]
    return(setu)
