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
##                2022-11-17 (QV) added 2D interpolation for a given sensor
##                2023-08-03 (QV) get lut url from ac.config
##                2024-02-26 (QV) added new lut support (with pressure)

def wvlut_interp(ths, thv, uwv=1.5, sensor=None, config='201710C',
                 par_id = 2, ## old LUT ttwava index
                 pressure = 1013, par = 'ttwava', ## new LUT, index computed
                 remote_base = None):
    import os, sys
    import scipy.interpolate
    import acolite as ac
    import numpy as np

    ## use URL from main config
    if remote_base is None: remote_base = '{}'.format(ac.config['lut_url'])

    ## input geometry dimensions
    dim = np.atleast_1d(ths).shape
    dim2 = np.atleast_1d(uwv).shape
    onedim = ((len(dim) == 1) & (dim[0] == 1)) & ((len(dim2) == 1) & (dim2[0] == 1))

    lut_path = '{}/LUT/WV'.format(ac.config['data_dir'])
    lut_id = 'WV_{}'.format(config)
    lutnc = '{}/{}.nc'.format(lut_path,lut_id)

    ## try downloading LUT from GitHub
    if (not os.path.isfile(lutnc)):
        remote_lut = '{}/WV/{}'.format(remote_base, os.path.basename(lutnc))
        try:
            print('Getting remote LUT {}'.format(remote_lut))
            ac.shared.download_file(remote_lut, lutnc)
            print('Testing LUT {}'.format(lutnc))
            lut, meta = ac.shared.lutnc_import(lutnc) # test LUT
        except:
            print('Could not download remote lut {} to {}'.format(remote_lut, lutnc))
            if os.path.exists(lutnc): os.remove(lutnc)

    ## import LUT
    if os.path.exists(lutnc):
        lut, meta = ac.shared.lutnc_import(lutnc)
    else:
        print('Could not open WV LUT {}'.format(lutnc))
        sys.exit(1)

    # find RSR
    if sensor is not None:
        rsrd = ac.shared.rsr_dict(sensor=sensor)
        if sensor in rsrd:
            rsr, rsr_bands = rsrd[sensor]['rsr'], rsrd[sensor]['rsr_bands']
        else:
            print('Sensor {} RSR not found'.format(sensor))
            sys.exit(1)

    if onedim:
        ## interpolate hyperspectral dataset
        ## old LUT
        if config=='201710C':
            rgi = scipy.interpolate.RegularGridInterpolator([meta['ths'], meta['thv'], meta['wv'], range(3),
                                                             meta['wave']],lut,
                                                             bounds_error=False, fill_value=None)
            iw = rgi((ths, thv, uwv, par_id, meta['wave']))
        ## new LUT
        else:
            ipd = {p:pi for pi, p in enumerate(meta['par'])}
            rgi = scipy.interpolate.RegularGridInterpolator([meta['pressure'], range(len(meta['par'])),
                                                             meta['wave'], meta['vza'], meta['sza'], meta['wv']],lut,
                                                                 bounds_error=False, fill_value=None)
            iw = rgi((pressure, ipd[par], meta['wave'], thv, ths, uwv))

        if sensor is None:
            ## return hyperspectral dataset for this geometry
            return(meta["wave"], iw)
        else:
            ## make band averaged values
            band_averaged = ac.shared.rsr_convolute_dict(meta['wave'], iw, rsr)
            return(band_averaged)
    else:
        ## make band averaged values
        if sensor == None:
            print('Multidimensional WV LUT interpolation not currently supported for hyperspectral.')
            return()

        ## convolution lut and make rgi
        band_averaged = {}
        for band in rsr_bands:
            ## old LUT
            if config=='201710C':
                blut = ac.shared.rsr_convolute_nd(lut, meta['wave'],rsr[band]['response'], rsr[band]['wave'], axis=4)
                rgi = scipy.interpolate.RegularGridInterpolator([meta['ths'], meta['thv'], meta['wv'], range(3)],blut,\
                                                                bounds_error=False, fill_value=None)
                iw = rgi((ths, thv, uwv, par_id))
            ## new LUT
            else:
                blut = ac.shared.rsr_convolute_nd(lut, meta['wave'],rsr[band]['response'], rsr[band]['wave'], axis=2)
                ipd = {p:pi for pi, p in enumerate(meta['par'])}
                rgi = scipy.interpolate.RegularGridInterpolator([meta['pressure'], range(len(meta['par'])),
                                                                 meta['vza'], meta['sza'], meta['wv']],blut,
                                                                 bounds_error=False, fill_value=None)
                iw = rgi((pressure, ipd[par], thv, ths, uwv))

            band_averaged[band] = iw
        return(band_averaged)
