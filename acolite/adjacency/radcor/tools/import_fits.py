## def import_fits
## import APSFS fit datasets for RAdCor
##
## written by Alexandre Castagna, UGent (R version)
##            Quinten Vanhellemont, RBINS (Python version)
##
## for the RAdCor project 2023-2024
##
## modifications: 2024-04-24 (QV) added as function, added sensor aliases dict, added PS_22
##                2025-02-01 (QV) added Pleiades 1A and 1B, WorldView2 and 3
##                2025-02-05 (QV) added MERIS

def import_fits(sensor, bands, aer_models = ['C', 'M']):
    import acolite as ac
    import json

    if ac.settings['run']['verbosity'] > 2: print('Loading APSF and SAF fits for {}'.format(sensor))

    ## determine bands to import
    bands_ = [b for b in bands if bands[b]['radcor_use_band']]

    ## fit directories
    psffit_dir = '{}/{}'.format(ac.config['data_dir'], 'Shared/RAdCor/psf/fit/')    # Aerosol Point Spread Function
    saffit_dir = '{}/{}'.format(ac.config['data_dir'], 'Shared/RAdCor/saf/fit/')    # Aerosol Spherical albedo
    raypsf_dir = '{}/{}'.format(ac.config['data_dir'], 'Shared/RAdCor/psf_ray/fit') # Rayleigh PSF
    raysaf_dir = '{}/{}'.format(ac.config['data_dir'], 'Shared/RAdCor/saf_ray/fit') # Rayleigh SAF

    ## sensor aliases to APSFS outputs
    aliases = {'PlanetScope_SD8': 'SD_00', 'PlanetScope_SD5': 'SD_00',
               'PlanetScope_0c': 'PS_0C0D', 'PlanetScope_0d05': 'PS_0C0D', 'PlanetScope_0d06': 'PS_0C0D',
               'PlanetScope_0e': 'PS_0E', 'PlanetScope_0f': 'PS_0F10', 'PlanetScope_22': 'PS_22',
               'PHR1A': 'Pld_PHRA', 'PHR1B': 'Pld_PHRB', 'WorldView2': 'WV2', 'WorldView3': 'WV3',
               'EN1_MERIS': 'EV_MERIS',}

    ## Load PSF data
    coefs_psf_ray = {}
    for fvar in ['fa_generic_R']:
        ffile = '{}/{}.json'.format(raypsf_dir, fvar)
        with open(ffile) as f:
            fd = json.load(f)
        coefs_psf_ray[fvar] = fd[fvar]

    coefs_psf_aer = {}
    for ai, am in enumerate(aer_models):
        aer_psf_f = {}
        for i, b in enumerate(bands_):
            if sensor in aliases:
                fvar = 'fa_{}_{}_band_{}'.format(aliases[sensor], am, bands[b]['fit_name'])
            else:
                fvar = 'fa_{}_{}_band_{}'.format(sensor, am, bands[b]['fit_name'])

            ffile = '{}/{}.json'.format(psffit_dir, fvar)
            with open(ffile) as f:
                fd = json.load(f)
            aer_psf_f[b] = fd[fvar]
        coefs_psf_aer[am] = aer_psf_f

    ## Load SAF data
    coefs_saf_ray = {}
    for fvar in ['fa_generic_R']:
        ffile = '{}/{}.json'.format(raysaf_dir, fvar)
        with open(ffile) as f:
            fd = json.load(f)
        coefs_saf_ray[fvar] = fd[fvar]

    coefs_saf_aer = {}
    for ai, am in enumerate(aer_models):
        aer_saf_f = {}
        for i, b in enumerate(bands_):
            if sensor in aliases:
                fvar = 'fa_{}_{}_band_{}'.format(aliases[sensor], am, bands[b]['fit_name'])
            else:
                fvar = 'fa_{}_{}_band_{}'.format(sensor, am, bands[b]['fit_name'])

            ffile = '{}/{}.json'.format(saffit_dir, fvar)
            with open(ffile) as f:
                fd = json.load(f)
            aer_saf_f[b] = fd[fvar]
        coefs_saf_aer[am] = aer_saf_f

    return(coefs_psf_ray, coefs_psf_aer, coefs_saf_ray, coefs_saf_aer)
