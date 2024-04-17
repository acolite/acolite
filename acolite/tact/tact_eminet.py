# def tact_eminet
## runs eminet for tact processing
##
## Newly trained nets based on https://doi.org/10.1016/j.isprsjprs.2020.06.007
## Now including L5, L7, L8, L9 and some extra filtering on ECOSTRESS spectra
##
## written by Quinten Vanhellemont, RBINS
## 2022-08-09
## modifications: 2022-08-11 (QV) don't attempt to run with missing bands
##                2022-10-23 (QV) added EMINET download
##                2024-04-17 (QV) use new gem NetCDF handling

def tact_eminet(gem, model_path = None,
                      water_fill = True, water_threshold = 0.0215,
                      fill = True, fill_dilate = False,
                      model_base = 'EMINET_{}_{}_64x4.h5', model_version = '20220809', netname = 'Net1', verbosity = 5):

    import os, json
    import acolite as ac
    from keras.models import load_model
    import numpy as np

    ## read gem file if NetCDF
    if type(gem) is str:
        gem = ac.gem.gem(gem)
    gemf = gem.file

    ## find L2R with rhos
    dn = os.path.dirname(gemf)
    bn = os.path.basename(gemf)
    ncf = '{}/{}'.format(dn, bn.replace('_L1R', '_L2R'))

    if not os.path.exists(ncf):
        print('ACOLITE L2R file required for EMINET {} not found.'.format(ncf))
        return(None)

    print('Opening {}'.format(ncf))
    gemi = ac.gem.gem(ncf)
    datasets = gemi.datasets
    sensors = { 'L5_TM':'L5_TM_B6', 'L7_ETM':'L7_ETM_B6',
                'L8_OLI':'L8_TIRS', 'L9_OLI':'L9_TIRS'}
    satsen = sensors[gemi.gatts['sensor']]

    ## select model
    if model_path is None:
        ## model and meta file names
        model_file = model_base.format(satsen, netname)
        meta_file = model_file[0:-3]+'_meta.json'
        ## local file
        model_dir = ac.config['data_dir'] + '/EMINET'
        model_path = '{}/{}/{}'.format(model_dir, model_version, model_file)
        meta_path = '{}/{}/{}'.format(model_dir, model_version, meta_file)
        ## remote file
        url_base = 'https://github.com/acolite/acolite_extra/raw/main/EMINET/'
        model_url = '{}/{}/{}'.format(url_base, model_version, model_file)
        meta_url = '{}/{}/{}'.format(url_base, model_version, meta_file)
        if not os.path.exists(model_path):
            print('Getting remote file {}'.format(model_url))
            ac.shared.download_file(model_url, model_path)
        if not os.path.exists(meta_path):
            print('Getting remote file {}'.format(meta_url))
            ac.shared.download_file(meta_url, meta_path)

    ## open model
    if verbosity > 2: print('Opening model file {}'.format(model_path))
    model = load_model(model_path)

    ## read metadata
    mfile = model_path.replace('.h5', '_meta.json')
    if verbosity > 2: print('Opening model metadata {}'.format(mfile))
    meta = json.load(open(mfile, 'r'))

    if water_fill:
        ## load water emissivity
        emissivity_file = '{}/{}/emissivity_{}.json'.format(ac.config['data_dir'], 'TACT', 'water')
        em = json.load(open(emissivity_file, 'r'))

    ## load RSR
    if satsen in ['L8_TIRS', 'L9_TIRS']:
        ibands = ['1','2','3','4','5','6','7']
        obands = ['10', '11']
        emsen = '{}'.format(satsen)

    if satsen in ['L5_TM_B6', 'L7_ETM_B6']:
        ibands = ['1','2','3','4','5','7']
        obands = ['6']
        emsen = '{}'.format(satsen[0:-3])
    rsr = ac.shared.rsr_dict(sensor = gemi.gatts['sensor'])[gemi.gatts['sensor']]
    datasets_needed = ['rhos_{}'.format(rsr['wave_name'][b]) for b in ibands]

    datasets_found = [ds for ds in datasets_needed if ds in datasets]
    if (datasets_found != datasets_needed):
        print('Could not load required datasets from {}'.format(ncf))
        print('Required: {}'.format(','.join(datasets_needed)))
        print('Found: {}'.format(','.join(datasets_found)))
        return(None)

    spect = None
    for ds in datasets_needed:
        data = gemi.data(ds)
        if len(np.where(np.isfinite(data))[0]) == 0: break
        if spect is None:
            data_dim = data.shape
            edge = np.isnan(data)
            spect = data.flatten()
        else:
            spect = np.vstack((spect, data.flatten()))
    gemi.close()
    gemi = None

    ## Maybe not needed if L1R is empty L2R will be missing
    if spect is None:
        print('Not running EMINET: Empty VSWIR data in {}'.format(ncf))
        return(None)

    ## transpose to n,7
    spect = spect.T

    ## use water emissivity defaults for pixels with low swir reflectance
    if water_fill:
        water_sub = np.where(spect[:, -2] < water_threshold)[0]
        print('Found {} "water" pixels'.format(len(water_sub)))

    ## reflectance to percent
    spect *= 100

    ## normalise to input distribution
    for ik in range(spect.shape[1]):
        spect[:,ik] = (spect[:,ik]-meta['xmean'][ik])/meta['xstd'][ik]

    ## run prediction
    if verbosity > 0: print('Running emissivity retrieval for {} ({} spectra)'.format(netname, spect.shape[0]))
    tir_pred = model.predict(spect)

    ## mask reflectance <0 and >100
    tir_pred[tir_pred<0] = np.nan
    tir_pred[tir_pred>100] = np.nan

    ## convert reflectance to emissivity 0--1
    em_pred = 1-tir_pred/100
    tir_pred = None

    result = {}
    for ki, k in enumerate(obands):
        result[k] = em_pred[:, ki]

        ## set water values
        if water_fill:
            result[k][water_sub] = em[emsen][k]

        ## reform to 2D
        result[k] = result[k].reshape(data_dim)

        ## fill holes with nearest data
        if fill:
            if fill_dilate:
                import scipy.ndimage
                struct = scipy.ndimage.generate_binary_structure(2, 2)
                tmp = scipy.ndimage.binary_dilation(np.isnan(result[k]), structure=struct, iterations=5)
                result[k][tmp] = np.nan
            ## fill nans
            result[k] = ac.shared.fillnan(result[k])

        ## mask edges
        result[k][edge] = np.nan
    return(result)
