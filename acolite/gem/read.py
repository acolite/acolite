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
        gem['gatts']['thermal_bands'] = ['6_vcid_1', '6_vcid_2']

    if 'projection_key' in gem['gatts']:
        gem['nc_projection'] = ac.shared.nc_read_projection(ncf)

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
