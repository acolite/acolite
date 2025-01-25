## read L1R.nc gem file
## written by Quinten Vanhellemont, RBINS
## 2021-02-28
## modifications: 2021-03-09 (QV) made reading data optional
##                2021-12-08 (QV) added nc_projection
##                2022-02-15 (QV) added L9/TIRS

def read(ncf, sub = None, skip_datasets = [], load_data=True):
    import os
    import numpy as np
    import acolite as ac

    if not os.path.exists(ncf): return()

    ## set up dict
    gem = {'data':{}, 'atts':{}}

    ## get datasets and attributes from NetCDF
    gem['datasets'] = ac.shared.nc_datasets(ncf)
    gem['gatts'] = ac.shared.nc_gatts(ncf)
    gem['gatts']['gemfile'] = ncf

    ## detect thermal sensor
    if gem['gatts']['sensor'] == 'L8_OLI':
        gem['gatts']['thermal_sensor'] = 'L8_TIRS'
        gem['gatts']['thermal_bands'] = ['10', '11']
    elif gem['gatts']['sensor'] == 'L9_OLI':
        gem['gatts']['thermal_sensor'] = 'L9_TIRS'
        gem['gatts']['thermal_bands'] = ['10', '11']
    elif gem['gatts']['sensor'] == 'L5_TM':
        gem['gatts']['thermal_sensor'] = 'L5_TM'
        gem['gatts']['thermal_bands'] = ['6']
    elif gem['gatts']['sensor'] == 'L7_ETM':
        gem['gatts']['thermal_sensor'] = 'L7_ETM'
        gem['gatts']['thermal_bands'] = ['6_vcid_1', '6_vcid_2', '6_VCID_1', '6_VCID_2']
    elif gem['gatts']['sensor'] == 'SUOMI-NPP_VIIRS':
        gem['gatts']['thermal_sensor'] = 'SUOMI-NPP_VIIRS_TIR'
        gem['gatts']['thermal_bands'] = ['I04', 'I05', 'M12', 'M13', 'M14', 'M15', 'M16']
    elif gem['gatts']['sensor'] == 'JPSS-1_VIIRS':
        gem['gatts']['thermal_sensor'] = 'JPSS-1_VIIRS_TIR'
        gem['gatts']['thermal_bands'] = ['I04', 'I05', 'M12', 'M13', 'M14', 'M15', 'M16']
    elif gem['gatts']['sensor'] == 'JPSS-2_VIIRS':
        gem['gatts']['thermal_sensor'] = 'JPSS-2_VIIRS_TIR'
        gem['gatts']['thermal_bands'] = ['I04', 'I05', 'M12', 'M13', 'M14', 'M15', 'M16']

    if 'projection_key' in gem['gatts']:
        gem['nc_projection'] = ac.shared.nc_projection_read(ncf)

    if load_data:
        ## read all datasets
        for ds in gem['datasets']:
            if ds in skip_datasets: continue
            if 'projection_key' in gem['gatts']:
                if ds in ['x', 'y', gem['gatts']['projection_key']]: continue
            d_, a_ = ac.shared.nc_data(ncf, ds, sub=sub, attributes=True)
            gem['data'][ds] = d_.data
            gem['data'][ds][d_.mask] = np.nan
            gem['atts'][ds] = a_

    return(gem)
