## o2lut_interp
## gives O2 transmittance for given sun and view zenith angles, sensor
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-10-30
## modifications:
##                2021-02-24 (QV) new interpolation, lut is determined here and read generically

def o2lut_interp(ths, thv, sensor=None, o2config='201810C', par_id = 2):
    import os
    import scipy.interpolate
    import acolite as ac

    lut_path = '{}/LUT/O2'.format(ac.config['data_dir'])
    lut_id = 'O2_{}'.format(o2config)
    lutnc = '{}/{}.nc'.format(lut_path,lut_id)
    lut, meta = ac.shared.lutnc_import(lutnc)

    ## interpolate hyperspectral dataset
    rgi = scipy.interpolate.RegularGridInterpolator([meta['ths'], meta['thv'], range(3), meta['wave']],lut,
                                                     bounds_error=False, fill_value=None)
    iw = rgi((ths, thv, par_id, meta['wave']))

    if sensor is None:
        ## return hyperspectral dataset for this geometry
        return(meta["wave"], iw)
    else:
        # find RSR
        rsr_file = ac.config['data_dir']+'/RSR/{}.txt'.format(sensor)
        rsr,bands = ac.shared.rsr_read(file=rsr_file)

        ## make band averaged values
        band_averaged = ac.shared.rsr_convolute_dict(meta['wave'], iw, rsr)
        return(band_averaged)
