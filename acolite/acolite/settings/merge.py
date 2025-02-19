## def merge
## merges run settings with sensor and additional settings
## returns merged "setu" settings
## written by Quinten Vanhellemont, RBINS
## 2025-02-11
## modifications:

def merge(sensor = None, settings = None):
    import acolite as ac

    ## get run settings
    ## should have been updated with user settings in acolite_run
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## get sensor specific defaults
    if sensor is not None:
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
    ## end set sensor specific defaults

    ## additional settings for run
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    return(setu)
