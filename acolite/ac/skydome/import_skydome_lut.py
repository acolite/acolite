## def import_skydome_lut
## function to import skydome lut
## returns lut (or lut dict if sensor specified) and lut metadata, with rgi
## written by Quinten Vanhellemont, RBINS
## 2025-10-30
## modifications:

def import_skydome_lut(lutnc = None, lut_base_skydome = 'ACOLITE-FFSS-202510-82W', sensor = None, override = False,
                       return_rgi = True, include_wave_dim = False,
                       par = 'rsky', get_remote = True, remote_base = None):
    import os
    import numpy as np
    import scipy.interpolate
    import acolite as ac

    rgi = None

    ## use URL from main config
    if remote_base is None: remote_base = '{}'.format(ac.config['lut_url'])
    if lut_base_skydome is None: lut_base_skydome = ac.settings['run']['lut_base_skydome']

    if lutnc is None:
        lutdir = '{}/{}/{}'.format(ac.config['lut_dir'], 'SKYDOME', lut_base_skydome)
        lutnc = '{}-{}.nc'.format(lutdir, par)
    else:
        if not os.path.exists(lutnc):
            print('File does not exist: {}'.format(lutnc))
            return
        lutdir = os.path.dirname(lutnc)

    if sensor is None:
        ## get remote generic LUT
        if (not os.path.isfile(lutnc)) & (get_remote):
            remote_lut = '{}/{}/{}'.format(remote_base, 'SKYDOME', os.path.basename(lutnc))
            try:
                ac.shared.download_file(remote_lut, lutnc)
            except:
                print('Could not download remote lut {} to {}'.format(remote_lut, lutnc))

        ## import generic LUT
        lut, meta = ac.shared.lutnc_import(lutnc)

        ## set up interpolator
        if return_rgi:
            ## last dimension of LUT is wavelength, can either be included in the interpolator or not
            if include_wave_dim:
                dim = [np.atleast_1d(meta[d]) for d in meta['dims']]
            else:
                dim = [np.atleast_1d(meta[d]) for d in meta['dims'][0:-1]]
            ## set up rgi
            rgi = scipy.interpolate.RegularGridInterpolator(dim, lut, bounds_error = False, fill_value = np.nan)
    else:
        slutdir = '{}/{}/{}/{}'.format(ac.config['lut_dir'], 'SKYDOME', sensor, lut_base_skydome)
        slutnc = '{}-{}_{}.nc'.format(slutdir, par, sensor)

        ## try downloading sensor LUT from GitHub
        if (not os.path.exists(slutnc)) & (get_remote):
            remote_lut = '{}/{}/{}/{}'.format(remote_base, 'SKYDOME', sensor, os.path.basename(slutnc))
            try:
                ac.shared.download_file(remote_lut, slutnc)
            except:
                print('Could not download remote lut {} to {}'.format(remote_lut, slutnc))
        rsrd = None

        ## create sensor LUT here
        if not os.path.exists(slutnc):
            print('Reading input LUT: {}'.format(lutnc))
            ## import generic LUT
            lut, meta = ac.shared.lutnc_import(lutnc)

            ## load spectral response
            rsrd = ac.shared.rsr_dict(sensor = sensor)

            ## resample to bands
            print('Resampling to sensor: {}'.format(sensor))
            lut_sensor = {}
            for band in rsrd[sensor]['rsr_bands']:
                lut_sensor[band] = ac.shared.rsr_convolute_nd(lut, meta['wavelength'],
                                                                   rsrd[sensor]['rsr'][band]['response'], rsrd[sensor]['rsr'][band]['wave'],
                                                                   axis = len(lut.shape)-1)

            ## write sensor LUT
            print('Writing output LUT: {}'.format(slutnc))
            metas = {k: meta[k] for k in meta if k not in ['dims', 'wavelength']}
            metas['dims'] = meta['dims'][0:-1]
            metas['sensor'] = sensor
            del lut, meta
            ac.shared.lutnc_write(slutnc, lut_sensor, metas)
            del lut_sensor, metas

        ## read sensor lut
        if os.path.exists(slutnc):
            lut, meta = ac.shared.lutnc_import(slutnc)

            if return_rgi:
                ## wavelength dim is not included for sensors
                dim = [np.atleast_1d(meta[d]) for d in meta['dims']]
                rgi = {}
                for band in lut:
                    rgi[band] = scipy.interpolate.RegularGridInterpolator(dim, lut[band], bounds_error = False, fill_value = np.nan)

    ## return metadata, lut, and rgi
    if return_rgi:
        return(meta, lut, rgi)

    ## return metadata and lut
    return(meta, lut)
