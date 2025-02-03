## def merged_lut
## separate function to merge atmosphere/interface LUT, saved to NetCDF if store = True
## written by Quinten Vanhellemont, RBINS
## 2024-04-25
## modifications:

def merged_lut(lut, lutint, pressure, sensor = None, get_remote = True, lut_par = None, store = True, override = False):
    import os
    import numpy as np
    import acolite as ac

    ## base lut
    lutid = '{}-{}mb'.format(lut, '{}'.format(pressure).zfill(4))
    lutdir = '{}/{}'.format(ac.config['lut_dir'], '-'.join(lutid.split('-')[0:3]))

    ## model
    model = lut.split('-MOD')[1] # lut[-1]
    mlut = lutid + lutint.replace('ACOLITE', '')
    mlut_nc = '{}/{}.nc'.format(lutdir, mlut)

    ## merged lut name
    if sensor is not None:
        mlut_nc = '{}/{}/{}_{}.nc'.format(lutdir, sensor, mlut, sensor)

    if (not os.path.exists(mlut_nc)) | (override):
        print('Generating merged LUT {}'.format(os.path.basename(mlut_nc)))
        lut_data, lut_meta = ac.aerlut.import_lut(lutid, lutdir, sensor = sensor, lut_par = None, get_remote = get_remote)
        rskyd = ac.aerlut.import_rsky_luts(sensor = sensor, models = [model], lutbase = lutint, get_remote = get_remote)

        rlut = rskyd[model]['lut']
        if 'wind' in rskyd[model]['meta']:
            rsky_winds  = rskyd[model]['meta']['wind']
        else:
            rsky_winds = np.atleast_1d([0,20])

        ## current pars
        ipd = {p:i for i,p in enumerate(lut_meta['par'])}

        ## sensor specific LUT
        if sensor is not None:
            dim_names = ['par', 'azi', 'thv', 'ths', 'wnd', 'tau']
            rsr_file = ac.config['data_dir']+'/RSR/'+sensor+'.txt'
            rsr, rsr_bands = ac.shared.rsr_read(file=rsr_file)
            lut_meta['bands'] = rsr_bands

            ## run through bands
            for band in rsr_bands:
                tlut = lut_data[band]
                tmp = rlut[band] * 1.0

                ## add wind axis
                if 'wind' not in rskyd[model]['meta']:
                    tmp = np.expand_dims(tmp, -2)
                    tmp = np.repeat(tmp, len(rsky_winds), axis=-2)
                tlut = np.repeat(tlut, len(rsky_winds), axis=-2)

                ## add to the LUT
                ## model rsky at surface
                ax = len(ipd)
                tlut = np.insert(tlut, (ax), tmp, axis=0) ## axis = 0 if pressure not in lut, first axis are params

                ## model rsky at toa
                ## (utott * dtott * rsky) / (1. - rsky * astot)
                ## first index parameters
                tmp = (tlut[ipd['utott'],:,:,:,:,:]*\
                        tlut[ipd['dtott'],:,:,:,:,:]*
                        tlut[ax,:,:,:,:,:]) /\
                         (1.-tlut[ax,:,:,:,:,:] *\
                         tlut[ipd['astot'],:,:,:,:,:])
                tlut = np.insert(tlut, (ax+1), tmp, axis=0)## axis = 0 if pressure not in lut, first axis are params

                ## add romix+rsky
                ## first index parameters
                tmp = tlut[ipd['romix'],:,:,:,:,:] + tlut[ax+1,:,:,:,:,:]
                tlut = np.insert(tlut, (ax+2), tmp, axis=0) ## axis = 0 if pressure not in lut, first axis are params

                ## replace in dict
                lut_data[band] = tlut

        ## generic LUT
        else:
            dim_names = ['par', 'wave', 'azi', 'thv', 'ths', 'wnd', 'tau']

            ## repeat 6sv lut for winds
            tlut = np.repeat(lut_data, len(rsky_winds), axis=-2)

            ## add to the LUT
            ## model rsky at surface
            ax = len(ipd)
            tlut = np.insert(tlut, (ax), rlut, axis=0) ## axis = 0 if pressure not in lut, first axis are params

            ## model rsky at toa
            ## (utott * dtott * rsky) / (1. - rsky * astot)
            tmp = (tlut[ipd['utott'],:,:,:,:,:]*\
                    tlut[ipd['dtott'],:,:,:,:,:]*
                    tlut[ax,:,:,:,:,:]) /\
                    (1.-tlut[ax,:,:,:,:,:] *\
                    tlut[ipd['astot'],:,:,:,:,:])
            tlut = np.insert(tlut, (ax+1), tmp, axis=0) ## axis = 0 if pressure not in lut, first axis are params

            ## add romix+rsky
            tmp = tlut[ipd['romix'],:,:,:,:,:] + tlut[ax+1,:,:,:,:,:]
            tlut = np.insert(tlut, (ax+2), tmp, axis=0) ## axis = 0 if pressure not in lut, first axis are params

            ## replace lut data
            lut_data = tlut

        del tmp, tlut

        ## add new pars
        lut_meta['par'] += ['rsky_s', 'rsky_t', 'romix+rsky_t']
        lut_meta['wnd'] = rsky_winds

        ## we'll store without compression for speed
        if store:
            #print('Writing {}'.format(os.path.basename(mlut_nc)))
            ac.shared.lutnc_write(mlut_nc, lut_data, lut_meta, dims = dim_names, compression = False)
            #ac.shared.lutnc_write(mlut_nc, lut_data, lut_meta, dims = dim_names, compression = True, complevel = 1)

    else:
        print('Reading {}'.format(os.path.basename(mlut_nc)))
        lut_data, lut_meta = ac.shared.lutnc_import(mlut_nc)

    ## subset to parameters
    if lut_par is not None:
        lut_sub_idx = []
        lut_sub_par = []
        for i, ik in enumerate(lut_par):
            for j, jk in enumerate(lut_meta['par']):
                if (ik == jk):
                    lut_sub_idx.append(j)
                    lut_sub_par.append(jk)
        if sensor is not None:
            for b in lut_data:
                lut_data[b] = lut_data[b][lut_sub_idx,:,:,:,:,:]
        else:
           lut_data = lut_data[lut_sub_idx,:,:,:,:,:]
        lut_meta['par'] = lut_sub_par

    return(lut_data, lut_meta)
