## def import_skydome_lut_reformat
## temporary function to reformat skydome to  LUT dimensions
## written by Quinten Vanhellemont, RBINS
## 2025-11-01
## modifications:

def import_skydome_lut_reformat(sensor = None, model = 1, lut_base_skydome = 'ACOLITE-FFSS-202511-82W'):
    import numpy as np
    import scipy.interpolate

    import acolite as ac

    ## fake wind dimension
    rsky_winds = np.atleast_1d(0)

    ## new raa and vza dims
    raa_new = np.array([ 0.,  10.,  20.,  40.,  60.,  80.,  90., 100., 120., 140., 160., 170., 180.])
    vza_new = np.array([ 0. ,  1.5,  4. ,  8. , 12. , 16. , 24. , 32. , 40. , 48. , 56. , 64. , 72. ])

    ## read skydome
    meta_rsky, lut_rsky = ac.ac.skydome.import_skydome_lut(sensor = sensor, par = 'rsky', lut_base_skydome = lut_base_skydome, return_rgi = False)

    ## compute fresnel
    muv = np.radians(meta_rsky['qza'])
    muv[muv == 0.0] = 0.1
    rhof = ac.ac.sky_refl(muv, n_w=1.34)

    ## reaarange and interpolate LUT
    if sensor is None:
        ## dimensions are ordered ['sza', 'qza', 'raa', 'aerosol', 'tau', 'wave']
        ## index model
        blut = lut_rsky[:,:, :, model, :, :]

        ## convert to fresnel reflected sky reflectance
        blut *= rhof[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis]

        ## interpolator to original LUT raa dimensions
        inta = scipy.interpolate.interp1d(np.abs(180 - meta_rsky['raa']), blut, axis=2) ## flip raa here to align skydome elements
        blut = inta(raa_new)

        ## interpolator to original LUT vza dimensions
        inta = scipy.interpolate.interp1d(meta_rsky['qza'], blut, axis=1)
        blut = inta(vza_new)

        ## reorder axis
        ## dimensions were sza, vza, raa, aot, wave
        ## dimensions have to be wave, raa, vza, sza, wind, aot
        blut = np.moveaxis(blut, [0, 1, 2, 3, 4], [3, 2, 1, 4, 0])

        ## add empty wind dim
        blut = np.expand_dims(blut, -2)
        lut_out = np.repeat(blut, len(rsky_winds), axis=-2)
        del blut
    else:
        lut_out = {}
        for b in lut_rsky:
            ## dimensions are ordered ['sza', 'qza', 'raa', 'aerosol', 'tau']
            ## index model
            blut = lut_rsky[b][:,:, :, model, :]

            ## convert to fresnel reflected sky reflectance
            blut *= rhof[np.newaxis, :, np.newaxis, np.newaxis]

            ## interpolator to original LUT raa dimensions
            inta = scipy.interpolate.interp1d(np.abs(180 - meta_rsky['raa']), blut, axis=2) ## flip raa here to align skydome elements
            blut = inta(raa_new)

            ## interpolator to original LUT vza dimensions
            inta = scipy.interpolate.interp1d(meta_rsky['qza'], blut, axis=1)
            blut = inta(vza_new)

            ## reorder axis
            ## dimensions were sza, vza, raa, aot
            ## dimensions have to be raa, vza, sza, wind, aot
            blut = np.moveaxis(blut, [0, 1, 2, 3], [2, 1, 0, 3])

            ## add empty wind dim
            blut = np.expand_dims(blut, -2)
            lut_out[b] = np.repeat(blut, len(rsky_winds), axis=-2)
            del blut

    meta_out = {}
    meta_out['raa'] = raa_new
    meta_out['vza'] = vza_new
    meta_out['sza'] = meta_rsky['sza']
    meta_out['wind'] = rsky_winds

    return(meta_out, lut_out)
