## def crat_nc
## computes Ruddick et al. 2001 CRAT chlorophyll for a NetCDF file
## https://doi.org/10.1364/AO.40.003575
##
## written by Quinten Vanhellemont, RBINS
## 2026-06-10
## modifications:

def crat_nc(ncf, parameter = 'rhow', datasets = None, crat_config = 'defaults'):
    import os
    import acolite as ac
    import numpy as np

    parameter_list = ['Rrs', 'rhow', 'rhos', 'rhosu']
    sensor_list = ['PACE_OCI', 'PRISMA', 'ENMAP_HSI', 'DESIS_HSI']

    if parameter not in parameter_list:
        print('Parameter {} not configured, use one of: '.format(parameter))
        print(', '.join(parameter_list))
        return

    ## read data with xarray
    data, gatts = ac.shared.xr.read_rho(ncf, parameter = parameter, datasets = datasets)

    ## check if we have any of the required bands
    if data[parameter].shape == (0,):
        print('Parameter {} not in {} '.format(parameter, ncf))
        del data, gatts
        return

    ## check if we have any of the required bands
    if gatts['sensor'] not in sensor_list:
        print('Sensor {} not supported for CRAT '.format(gatts['sensor']))
        del data, gatts
        return

    ## get inputfile info from gatts
    bn = os.path.basename(gatts['inputfile'])
    version = bn[bn.find('V'):].split('.')[0:-1]
    sp = bn.split('.')

    ## compute crat
    if parameter == 'Rrs':
        c = ac.parameters.chl_crat.crat(data['wavelength'], data['Rrs'] * np.pi, crat_config = crat_config)
    elif parameter in ['rhow', 'rhos', 'rhosu']:
        c = ac.parameters.chl_crat.crat(data['wavelength'], data[parameter] * np.pi, crat_config = crat_config)

    if datasets is None:
        del data, gatts
        return(c)
    return(c, data, gatts)
