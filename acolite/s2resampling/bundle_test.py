## def bundle_test
## tests netcdf file for S2Resampling cues
##
## written by Quinten Vanhellemont, RBINS
## 2023-05-02
## modifications: 2023-05-05 (QV) update for HROC mosaics
##                2026-09-22 (QV) added S2C support

def bundle_test(file):
    import os
    import acolite as ac

    ## read gatts
    gatts = ac.shared.nc_gatts(file)

    ## identify sensor
    sensor = None
    if ('platform' in gatts) & ('type' in gatts): ## HROC mosaic
        if gatts['type'] == 'MSIL1C':
            sensor = '{}_MSI'.format(gatts['platform'])
    if sensor is None:
        dn = os.path.dirname(file)
        bn = os.path.basename(file)
        if bn[0:3] == 'S2A':
            sensor = 'S2A_MSI'
        elif bn[0:3] == 'S2B':
            sensor = 'S2B_MSI'
        elif bn[0:3] == 'S2C':
            sensor = 'S2C_MSI'
        else:
            return

    ## start date
    if 'start_date' not in gatts:
        return

    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## get sensor specific defaults
    setd = ac.acolite.settings.parse(sensor)
    ## set sensor default if user has not specified the setting
    for k in setd:
        if k not in ac.settings['user']: setu[k] = setd[k]
    ## end set sensor specific defaults

    ## sensor rsrd
    if setu['rsr_version'] is None:
        sensor_lut = '{}'.format(sensor)
    else:
        sensor_lut = '{}_{}'.format(sensor, setu['rsr_version'])
    rsrd = ac.shared.rsr_dict(sensor_lut)[sensor_lut]

    ## check datasets
    datasets = ac.shared.nc_datasets(file)
    required_datasets = ['lat', 'lon',  'view_zenith_mean', 'view_azimuth_mean', 'sun_zenith', 'sun_azimuth']
    for b in rsrd['rsr_bands']:
        ds = 'B{}'.format(b)
        required_datasets.append(ds)
    s2r_file = all([ds in datasets for ds in required_datasets])

    if s2r_file: return(sensor, gatts, datasets)
