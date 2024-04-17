# def ndvi_emissivity
## compute emissivity from NDVI for tact processing
##
## See e.g. Sobrino et al 2008, Skokovic et al 2014
## https://doi.org/10.1016/j.isprsjprs.2020.06.007 and references therein
##
## written by Quinten Vanhellemont, RBINS
## 2022-08-10
## modifications: 2024-04-17 (QV) removed keras import, use new gem NetCDF handling

def ndvi_emissivity(gem, ndvi_toa=True):

    import os, json
    import acolite as ac
    import numpy as np

    ## read gem file if NetCDF
    if type(gem) is str:
        gem = ac.gem.gem(gem)
    gemf = gem.file

    ## find L2R with rhos
    if ndvi_toa is False:
        ds_base = 'rhos_{}'
        dn = os.path.dirname(gemf)
        bn = os.path.basename(gemf)
        ncf = '{}/{}'.format(dn, bn.replace('_L1R', '_L2R'))
        if not os.path.exists(ncf):
            print('ACOLITE L2R file required for surface level NDVI {} not found.')
            return(None)
        open_l2r = True
    else:
        ds_base = 'rhot_{}'
        ncf = '{}'.format(gemf)
        open_l2r = False

    sensors = { 'L5_TM':'L5_TM_B6','L7_ETM':'L7_ETM_B6',
                'L8_OLI':'L8_TIRS', 'L9_OLI':'L9_TIRS'}

    ## open reflectance file
    if open_l2r:
        print('Opening {}'.format(ncf))
        gemi = ac.gem.gem(ncf)
        satsen = sensors[gemi.gatts['sensor']]
        datasets = gemi.datasets
    else:
        print('Using {}'.format(gemf))
        satsen = sensors[gem.gatts['sensor']]
        datasets = gem.datasets

    ## available wavelengths
    ds_waves = [ds.split('_')[1] for ds in datasets if ds_base[0:4]==ds[0:4]]

    ## find datasets for NDVI
    ## wavelengths and max wavelength difference
    ndvi_diff = [20, 40]
    req_waves = [660,865]
    datasets_needed = []
    for i, reqw in enumerate(req_waves):
        widx,selwave = ac.shared.closest_idx(ds_waves, reqw)
        if abs(float(selwave)-float(reqw)) > ndvi_diff[i]: continue
        datasets_needed.append(ds_base.format(selwave))
    datasets_found = [ds for ds in datasets_needed if ds in datasets]
    if (datasets_found != datasets_needed):
        print('Could not load required datasets from {}'.format(ncf))
        print('Required: {}'.format(','.join(datasets_needed)))
        print('Found: {}'.format(','.join(datasets_found)))
        return(None)

    ## compute ndvi
    if open_l2r:
        red = gemi.data(datasets_found[0])
        nir = gemi.data(datasets_found[1])
        gemi.close()
        gemi = None
    else:
        red = gem.data(datasets_found[0])
        nir = gem.data(datasets_found[1])

    ndvi = (nir-red)/(nir+red)

    if len(np.where(np.isfinite(ndvi))[0]) == 0:
        print('Not computing NDVI based emissivity: Empty VSWIR data in {}'.format(ncf))
        return(None)

    ndvi[np.abs(ndvi)>1]=np.nan
    nir = None

    if satsen in ['L8_TIRS', 'L9_TIRS']:
        ## Skokovic 2014
        ndvi_range = (0.18, 0.85)
        em_soil = {'10': 0.971, '11': 0.977}
        em_veg = {'10': 0.987, '11': 0.989}
        em_water = {'10': 0.991, '11': 0.986}
        em_ice = {'10': 0.986, '11': 0.959}
        fvc0 = {'10':[0.979, -0.046], '11':[0.982, - 0.027]}
    elif satsen in ['L5_TM_B6', 'L7_ETM_B6']:
        ## Sobrino 2008
        ndvi_range = (0.2, 0.5)
        em_veg = {'6': [0.004, 0.986]}
        fvc0 = {'6': [0.979, -0.035]}

    ## compute fvc
    fvc = ((ndvi-ndvi_range[0]) / (ndvi_range[1]-ndvi_range[0]))**2
    fvc[ndvi<ndvi_range[0]] = 0
    fvc[ndvi>ndvi_range[1]] = 1.

    em = {}
    em['fvc'] = fvc
    em['ndvi'] = ndvi

    for b in em_soil:
        if b == '6':
            em[b] = em_veg[b][0] * fvc + em_veg[b][1]
        else:
            em[b] = em_soil[b] * (1-fvc) + em_veg[b] * fvc
        em[b][fvc == 0] = fvc0[b][0] + fvc0[b][1] * red[fvc == 0]

    fvc = None
    ndvi = None
    red = None
    return(em)
