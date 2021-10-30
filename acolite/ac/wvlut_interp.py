## wvlut_interp
## gives WV transmittance for given sun and view zenith angles, wv and sensor
## uwv = water vapour in g/cm2
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-10-17
## modifications: 2017-10-18 (QV) added default uwv
##                2017-11-28 (QV) moved PP data directory
##                2018-01-31 (QV) fixed return when no sensor is given
##                2018-07-18 (QV) changed acolite import name
##                2021-02-24 (QV) new interpolation, lut is determined here and read generically

def wvlut_interp(ths, thv, uwv=1.5, sensor=None, config='201710C', par_id = 2,
                  remote_base = 'https://raw.githubusercontent.com/acolite/acolite_luts/main'):
    import os, sys
    import scipy.interpolate
    import acolite as ac

    lut_path = '{}/LUT/WV'.format(ac.config['data_dir'])
    lut_id = 'WV_{}'.format(config)
    lutnc = '{}/{}.nc'.format(lut_path,lut_id)

    ## try downloading LUT from GitHub
    if (not os.path.isfile(lutnc)):
        remote_lut = '{}/WV/{}'.format(remote_base, os.path.basename(lutnc))
        try:
            ac.shared.download_file(remote_lut, lutnc)
        except:
            print('Could not download remote lut {} to {}'.format(remote_lut, lutnc))

    ## import LUT
    if os.path.exists(lutnc):
        lut, meta = ac.shared.lutnc_import(lutnc)
    else:
        print('Could not open WV LUT {}'.format(lutnc))
        sys.exit(1)

    ## interpolate hyperspectral dataset
    rgi = scipy.interpolate.RegularGridInterpolator([meta['ths'], meta['thv'], meta['wv'], range(3), meta['wave']],lut,
                                                     bounds_error=False, fill_value=None)
    iw = rgi((ths, thv, uwv, par_id, meta['wave']))

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
