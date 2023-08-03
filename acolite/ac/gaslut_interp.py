## gaslut_interp
## returns gas transmittances (-H2O -O3) for given sun and view zenith angles, sensor
##
## written by Quinten Vanhellemont, RBINS
## 2021-06-22
## modifications: 2022-11-17 (QV) added 2D interpolation for a given sensor, removed waves keyword
##                2023-08-03 (QV) get lut url from ac.config

def gaslut_interp(sza, vza, pressure = 1013, sensor = None,
                  lutconfig = '202106F', pars = ['ttdica','ttoxyg','ttniox','ttmeth'],
                  remote_base = None):
    import os, sys
    import acolite as ac
    from netCDF4 import Dataset
    import scipy.interpolate
    import numpy as np

    ## use URL from main config
    if remote_base is None: remote_base = '{}'.format(ac.config['lut_url'])

    ## input geometry dimensions
    dim = np.atleast_1d(sza).shape
    onedim = ((len(dim) == 1) & (dim[0] == 1))

    ## identify LUT file
    lut_path = '{}/Gas'.format(ac.config['lut_dir'])
    lut_id = 'Gas_{}'.format(lutconfig)
    lutnc = '{}/{}.nc'.format(lut_path,lut_id)

    ## try downloading LUT from GitHub
    if (not os.path.isfile(lutnc)):
        remote_lut = '{}/Gas/{}'.format(remote_base, os.path.basename(lutnc))
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
        print('Could not open gas LUT {}'.format(lutnc))
        sys.exit(1)

    # find RSR
    if sensor is not None:
        rsrd = ac.shared.rsr_dict(sensor=sensor)
        if sensor in rsrd:
            rsr, rsr_bands = rsrd[sensor]['rsr'], rsrd[sensor]['rsr_bands']
        else:
            print('Sensor {} RSR not found'.format(sensor))
            sys.exit(1)

    ## LUT parameter indices
    ipd = {p:pi for pi, p in enumerate(meta['par'])}

    if onedim:
        ## set up interpolator
        rgi = scipy.interpolate.RegularGridInterpolator([meta['pressure'], range(len(meta['par'])),
                                                         meta['wave'], meta['vza'], meta['sza']],lut,
                                                                 bounds_error=False, fill_value=None)

        ## interpolate to vza, sza, pressure
        tg = {}
        for par in pars:
            tg[par] = rgi((pressure, ipd[par], meta['wave'], vza, sza))
        tg['wave'] = meta['wave']

        ## resample to sensor
        if sensor is not None:
            for par in tg:
                tg[par] = ac.shared.rsr_convolute_dict(meta['wave'], tg[par], rsr)
    else:
        ## make band averaged values
        if sensor == None:
            print('Multidimensional gas LUT interpolation not currently supported for hyperspectral.')
            return()

        ## convolution lut and make rgi
        tg = {}
        for band in rsr_bands:
            blut = ac.shared.rsr_convolute_nd(lut, meta['wave'],rsr[band]['response'], rsr[band]['wave'], axis=2)
            rgi = scipy.interpolate.RegularGridInterpolator([meta['pressure'], range(len(meta['par'])),
                                                             meta['vza'], meta['sza']],blut,
                                                             bounds_error=False, fill_value=None)
            ## run through parameters
            for par in ipd:
                if par not in tg: tg[par] = {}
                iw = rgi((pressure, ipd[par], vza, sza))
                tg[par][band] = iw

    return(tg)
