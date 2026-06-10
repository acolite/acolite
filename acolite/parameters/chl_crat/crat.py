## def crat
## computes Ruddick et al. 2001 CRAT chlorophyll
## https://doi.org/10.1364/AO.40.003575
##
## if rhow 1D or 2D, rhow is converted internally to 3D [wavelength, x, y]
##
## written by Quinten Vanhellemont, RBINS
## 2026-06-03
## modifications: 2026-06-10 (QV) added to ac.parameters.crat, changed settings to toml

def crat(wave, rhow, crat_config = 'defaults'):
    import acolite as ac
    import numpy as np
    import tomllib

    in_shape = rhow.shape
    ## add empty axes for 1 or 2D arrays
    if len(in_shape) != 3:
        new_axis = [1+i for i in range(3-len(in_shape))]
        rhow = np.expand_dims(rhow, axis = new_axis)

    ## read crat config
    crat_file = ac.config['data_dir'] + "/Shared/algorithms/chl_crat/crat.toml"
    crat_cfg = tomllib.load(open(crat_file, 'rb'))[crat_config]

    ## smooth input spectrum
    if crat_cfg['smooth']['iterations'] > 0:
        rhow = ac.shared.array.convolve(rhow, window = crat_cfg['smooth']['window'],
                                              iterations = crat_cfg['smooth']['iterations'])

    ## get water absorption
    awTS = ac.shared.wopp.aw_ts(crat_cfg['temperature'], crat_cfg['salinity'])

    ## compute reflectance at reference wavelength
    lw1 = ac.shared.array.linear_weight(wave, crat_cfg['wave']['ref'])
    rhow_ref = rhow[lw1[0], :, :] * lw1[1] + rhow[lw1[2], :, :] * lw1[3]

    ## find wavelengths where rhow == rhow_ref
    wsub = np.where((wave >= crat_cfg['wave']['min']) & (wave <= crat_cfg['wave']['max']))
    diff = rhow[wsub[0], :, :] - rhow_ref
    idx = np.argsort(np.abs(diff), axis = 0)[0, :, :]

    ## get wavelength indices
    s1 = wsub[0][idx]

    ## make sure first wavelength index >= 0
    s1i = s1-1
    s1i[s1i<0] = 0

    ## make sure last wavelength index <= len(wave)-1
    s1j = s1 + 1
    s1j[s1j > len(wave)-1] = len(wave)-1

    ## get bounding wavelengths for 2nd point
    wave_2_ = wave[s1]
    wave_2i = wave[s1i]
    wave_2j = wave[s1j]

    ## get bounding rhow for 2nd point
    rhow_2 = np.take_along_axis(rhow, np.expand_dims(s1, axis = 0), axis = 0)[0, :, :]
    rhow_2i = np.take_along_axis(rhow, np.expand_dims(s1i, axis = 0), axis = 0)[0, :, :]
    rhow_2j = np.take_along_axis(rhow, np.expand_dims(s1j, axis = 0), axis = 0)[0, :, :]

    ## interpolate from index before if rhow > rhow_ref
    out_2i = np.where(rhow_2 > rhow_ref)

    ## interpolate to index after if rhow < rhow_ref
    out_2j = np.where(rhow_2 < rhow_ref)

    ## create wave2 dataset
    wave_2 = wave_2_ * 1.0

    ## weighting for 2i
    vi = (rhow_ref[out_2i] - rhow_2i[out_2i]) / (rhow_2[out_2i] - rhow_2i[out_2i])
    ## these weights should not be >1 or <0 but this happens when the spectrum is not monotonous
    #vi[vi>1] = 1
    #vi[vi<0] = 0
    wave_2[out_2i] = wave_2_[out_2i] * (vi) + wave_2i[out_2i] * (1 - vi)

    ## weighting for 2j
    vj = (rhow_ref[out_2j] - rhow_2j[out_2j]) / (rhow_2[out_2j] - rhow_2j[out_2j])
    ## these weights should not be >1 or <0 but this happens when the spectrum is not monotonous
    #vj[vj>1] = 1
    #vj[vj<0] = 0
    wave_2[out_2j] = wave_2_[out_2j] * (vj) + wave_2j[out_2j] * (1 - vj)

    ## applicability mask
    tmp = (rhow_ref < rhow_2i) & (rhow_ref >= rhow_2j)
    wave_2[np.where(tmp == False)] = np.nan

    ## resample aw and compute a difference
    aw_1 = np.interp(crat_cfg['wave']['ref'], awTS['wave'], awTS['awTS'])
    aw_2 = np.interp(wave_2, awTS['wave'], awTS['awTS'])
    a = (aw_2 - aw_1)
    c = a / crat_cfg['aphy']

    ## NN version
    #aw_2nn = np.interp(wave_2_, awTS['wave'], awTS['awTS'])
    #ann = (aw_2nn-aw_1)
    #cnn = ann/ aphy
    #cnn[np.where(tmp == False)] = np.nan

    ## reshape output
    if len(in_shape) == 2: c = c.reshape(in_shape[1])

    ## reshape input
    if len(in_shape) != 3: rhow = rhow.reshape(in_shape)

    return(c)
