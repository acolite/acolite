## gaslut_interp
## returns gas transmittances (-H2O -O3) for given sun and view zenith angles, sensor
##
## written by Quinten Vanhellemont, RBINS
## 2021-06-22
## modifications:
##
def gaslut_interp(sza, vza, pressure = 1013,
                  sensor = None, waves = None,
                  lutconfig = '202106F', pars = ['ttdica','ttoxyg','ttniox','ttmeth'],
                  remote_base = 'https://raw.githubusercontent.com/acolite/acolite_luts/main'):
    import os, sys
    import acolite as ac
    from netCDF4 import Dataset
    import scipy.interpolate
    import numpy as np

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
        print('Could not open WV LUT {}'.format(lutnc))
        sys.exit(1)

    ## set up interpolator
    ipd = {p:pi for pi, p in enumerate(meta['par'])}
    rgi = scipy.interpolate.RegularGridInterpolator([meta['pressure'], range(len(meta['par'])),
                                                     meta['wave'], meta['vza'], meta['sza']],lut,
                                                             bounds_error=False, fill_value=None)

    ## interpolate to vza, sza, pressure
    tg = {}
    for par in pars:
        tg[par] = rgi((pressure, ipd[par], meta['wave'], vza, sza))
    tg['wave'] = meta['wave']

    ## interpolate to given wavelengths
    if waves is not None:
        for par in pars:
            tg[par] = np.interp(waves, meta['wave'], tg[par])

    ## resample to sensor
    if sensor is not None:
        # find RSR
        rsr_file = ac.config['data_dir']+'/RSR/{}.txt'.format(sensor)
        rsr,bands = ac.shared.rsr_read(file=rsr_file)
        ## make band averaged values
        for par in tg:
            tg[par] = ac.shared.rsr_convolute_dict(meta['wave'], tg[par], rsr)

    return(tg)
