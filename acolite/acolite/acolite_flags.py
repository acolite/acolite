## separate flagging function for ACOLITE NetCDF files
## flagging steps copied from acolite_l2w
##
## written by Quinten Vanhellemont, RBINS
## 2024-01-26
## modifications: 2024-05-21 (QV) skip negative rhos masking if rhos datasets are not present
##                2024-01-31 (QV) use first dataset to determine dimensions
##                2025-02-04 (QV) updated settings parsing

def acolite_flags(gem, create_flags_dataset=True, write_flags_dataset=False, return_flags_dataset=True):
    import acolite as ac
    import numpy as np
    import scipy.ndimage

    ## read gem file if NetCDF
    if type(gem) is str:
        gemf = '{}'.format(gem)
        gem = ac.gem.gem(gem)
    gemf = gem.file

    ## combine default and user defined settings
    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}
    ## get sensor specific defaults
    setd = ac.acolite.settings.parse(gem.gatts['sensor'])
    ## set sensor default if user has not specified the setting
    for k in setd:
        if k not in ac.settings['user']: setu[k] = setd[k]
    ## end set sensor specific defaults

    if setu['verbosity'] > 0: print('Running ACOLITE flagging function for {}.'.format(gemf))

    ## get rhot and rhos wavelengths
    rhot_ds = [ds for ds in gem.datasets if ds.startswith('rhot_')]
    rhot_waves = [int(ds.split('_')[-1]) for ds in rhot_ds]
    if len(rhot_waves) == 0: print('{} is probably not an ACOLITE L2R file: {} rhot datasets.'.format(gemf, len(rhot_ds)))

    rhos_ds = [ds for ds in gem.datasets if ds.startswith('rhos_')]
    rhos_waves = [int(ds.split('_')[-1]) for ds in rhos_ds]
    if len(rhos_waves) == 0: print('{} is probably not an ACOLITE L2R file: {} rhos datasets.'.format(gemf, len(rhos_ds)))

    ## source or output flags name
    ## only used when writing here
    flags_name = 'l2_flags'
    #flags_name = '{}_flags'.format(gem.gatts['acolite_file_type'][0:2].lower())

    ## create l2_flags dataset
    flags_att = {}
    if (flags_name not in gem.datasets) | (create_flags_dataset):
        dimensions = None
        for ds in gem.datasets:
            if ds in ['transverse_mercator', 'x', 'y',]: continue
            dimensions = gem.data(ds).shape
            break
        flags = np.zeros(dimensions,np.int32)
    else:
        flags = gem.data(flags_name)

    ## compute flags
    ####
    ## non water/swir threshold
    if setu['verbosity'] > 3: print('Computing non water threshold mask.')
    cidx,cwave = ac.shared.closest_idx(rhot_waves, setu['l2w_mask_wave'])
    ## use M bands for masking
    if ('VIIRS' in gem.gatts['sensor']) & (setu['viirs_mask_mband']):
        rhot_waves_m = [int(ds.split('_')[-1]) for ds in rhot_ds if 'M' in ds]
        cidx,cwave = ac.shared.closest_idx(rhot_waves_m, setu['l2w_mask_wave'])
    cur_par = 'rhot_{}'.format(cwave)
    cur_par = [ds for ds in rhot_ds if ('{:.0f}'.format(cwave) in ds)][0]
    if setu['verbosity'] > 3: print('Computing non water threshold mask from {} > {}.'.format(cur_par, setu['l2w_mask_threshold']))
    cur_data = gem.data(cur_par)
    if setu['l2w_mask_smooth']:
        cur_data = ac.shared.fillnan(cur_data)
        cur_data = scipy.ndimage.gaussian_filter(cur_data, setu['l2w_mask_smooth_sigma'], mode='reflect')
    cur_mask = cur_data > setu['l2w_mask_threshold']
    cur_data = None
    flags = cur_mask.astype(np.int32)*(2**setu['flag_exponent_swir'])
    cur_mask = None
    ## end non water/swir threshold
    ####

    ####
    ## cirrus masking
    if setu['verbosity'] > 3: print('Computing cirrus mask.')
    cidx,cwave = ac.shared.closest_idx(rhot_waves, setu['l2w_mask_cirrus_wave'])
    if np.abs(cwave - setu['l2w_mask_cirrus_wave']) < 5:
        cur_par = 'rhot_{}'.format(cwave)
        cur_par = [ds for ds in rhot_ds if ('{:.0f}'.format(cwave) in ds)][0]

        if setu['verbosity'] > 3: print('Computing cirrus mask from {} > {}.'.format(cur_par, setu['l2w_mask_cirrus_threshold']))
        cur_data = gem.data(cur_par)
        if setu['l2w_mask_smooth']:
            cur_data = ac.shared.fillnan(cur_data)
            cur_data = scipy.ndimage.gaussian_filter(cur_data, setu['l2w_mask_smooth_sigma'], mode='reflect')
        cirrus_mask = cur_data > setu['l2w_mask_cirrus_threshold']
        cirrus = None
        flags = (flags) | (cirrus_mask.astype(np.int32)*(2**setu['flag_exponent_cirrus']))
        cirrus_mask = None
    else:
        if setu['verbosity'] > 2: print('No suitable band found for cirrus masking.')
    ## end cirrus masking
    ####

    ####
    ## TOA out of limit
    if setu['verbosity'] > 3: print('Computing TOA limit mask.')
    toa_mask = None
    outmask = None
    for ci, cur_par in enumerate(rhot_ds):
        if rhot_waves[ci]<setu['l2w_mask_high_toa_wave_range'][0]: continue
        if rhot_waves[ci]>setu['l2w_mask_high_toa_wave_range'][1]: continue
        if setu['verbosity'] > 3: print('Computing TOA limit mask from {} > {}.'.format(cur_par, setu['l2w_mask_high_toa_threshold']))
        cwave = rhot_waves[ci]
        cur_par = [ds for ds in rhot_ds if ('{:.0f}'.format(cwave) in ds)][0]
        cur_data = gem.data(cur_par)
        if outmask is None: outmask = np.zeros(cur_data.shape).astype(bool)
        outmask = (outmask) | (np.isnan(cur_data))
        if setu['l2w_mask_smooth']:
            cur_data = ac.shared.fillnan(cur_data)
            cur_data = scipy.ndimage.gaussian_filter(cur_data, setu['l2w_mask_smooth_sigma'], mode='reflect')
        if toa_mask is None: toa_mask = np.zeros(cur_data.shape).astype(bool)
        toa_mask = (toa_mask) | (cur_data > setu['l2w_mask_high_toa_threshold'])
    flags = (flags) | (toa_mask.astype(np.int32)*(2**setu['flag_exponent_toa']))
    toa_mask = None
    flags = (flags) | (outmask.astype(np.int32)*(2**setu['flag_exponent_outofscene']))
    outmask = None
    ## end TOA out of limit
    ####

    ####
    ## negative rhos
    if len(rhos_waves) == 0:
        if setu['verbosity'] > 3: print('Not computing negative reflectance mask as there are no rhos datasets.')
    else:
        if setu['verbosity'] > 3: print('Computing negative reflectance mask.')
        neg_mask = None
        for ci, cur_par in enumerate(rhos_ds):
            if rhos_waves[ci]<setu['l2w_mask_negative_wave_range'][0]: continue
            if rhos_waves[ci]>setu['l2w_mask_negative_wave_range'][1]: continue
            if setu['verbosity'] > 3: print('Computing negative reflectance mask from {}.'.format(cur_par))
            cwave = rhos_waves[ci]
            cur_par = [ds for ds in rhos_ds if ('{:.0f}'.format(cwave) in ds)][0]
            cur_data = gem.data(cur_par)
            #if setu['l2w_mask_smooth']: cur_data = scipy.ndimage.gaussian_filter(cur_data, setu['l2w_mask_smooth_sigma'])
            if neg_mask is None: neg_mask = np.zeros(cur_data.shape).astype(bool)
            neg_mask = (neg_mask) | (cur_data < 0)
        flags = (flags) | (neg_mask.astype(np.int32)*(2**setu['flag_exponent_negative']))
        neg_mask = None
    ## end negative rhos
    ####

    ####
    ## mixed pixels mask
    if ('VIIRS' in gem.gatts['sensor']) & (setu['viirs_mask_immixed']):
        if setu['verbosity'] > 3: print('Finding mixed pixels using VIIRS I and M bands.')
        mix_mask = None
        if type(setu['viirs_mask_immixed_bands']) is not list:
            setu['viirs_mask_immixed_bands'] = [setu['viirs_mask_immixed_bands']]
        for imc in setu['viirs_mask_immixed_bands']:
            ib, mb = imc.split('/')
            ds0 = [ds for ds in rhot_ds if ib in ds][0]
            ds1 = [ds for ds in rhot_ds if mb in ds][0]
            cur_data0 = gem.data(ds0)
            cur_data1 = gem.data(ds1)
            if mix_mask is None: mix_mask = np.zeros(cur_data0.shape).astype(bool)
            if setu['viirs_mask_immixed_rat']:
                mix_mask = (mix_mask) | (np.abs(1-(cur_data0/cur_data1)) > setu['viirs_mask_immixed_maxrat'])
            if setu['viirs_mask_immixed_dif']:
                mix_mask = (mix_mask) | (np.abs(cur_data0-cur_data1) > setu['viirs_mask_immixed_maxdif'])
        flags = (flags) | (mix_mask.astype(np.int32)*(2**setu['flag_exponent_mixed']))
        mix_mask = None
    ## end mixed pixels mask
    ####

    ####
    ## dem shadow mask
    if setu['dem_shadow_mask']:
        if setu['verbosity'] > 3:
            print('Computing DEM shadow mask.')
            for k in setu:
                if 'dem_shadow' in k: print(k, setu[k])
        ## add gem version of dem_shadow_mask_nc?
        shade = ac.masking.dem_shadow_mask_nc(gemf)
        flags += shade.astype(np.int32)*(2**setu['flag_exponent_dem_shadow'])
    ## end dem shadow mask
    ####

    ## write flags dataset to netcdf
    if (write_flags_dataset): gem.write(flags_name, flags, ds_att = flags_att)

    ## return flags dataset
    if (return_flags_dataset): return(flags)

    return(gem)
