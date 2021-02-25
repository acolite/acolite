## read all luts and set up rgi
## QV 2020-01-15
## Last modifications: 2020-07-14 (QV) added sensor option
##                     2020-07-14 (QV) changed sensor option to a dict per band
##                     2020-07-22 (QV) added the rsky lut to the lut and rgi - need to update rsky azimuths
##                     2020-07-25 (QV) added updates to rsky lut
##                     2021-01-19 (QV) update to new LUTs
##                     2021-02-01 (QV) added wind speed for rsky lut, removed temp fixes for old luts
##                     2021-02-03 (QV) added down*up total transmittances

def import_luts(pressures = [500, 1013, 1100],
                base_luts = ['ACOLITE-LUT-202101-MOD1', 'ACOLITE-LUT-202101-MOD2'],
                rsky_lut = 'ACOLITE-RSKY-202101-75W-2ms',
                rsky_base = 'ACOLITE-RSKY-202101-75W-{}ms', rsky_winds = [1,2,5,10], add_rsky_winds = False,
                sensor=None, add_rsky=False, add_dutott = True):
    import scipy.interpolate
    import numpy as np
    import acolite as ac

    lut_dict = {}
    ## run through luts
    for lut in base_luts:
        ## run through pressures
        for ip, pr in enumerate(pressures):
            lutid = '{}-{}mb'.format(lut, '{}'.format(pr).zfill(4))
            lutdir = '{}/'.format(ac.config['lut_dir'])
            if sensor is None:
                ## LUT with 18 monochromatic wavelengths (0.39-2.4)
                lut_data, lut_meta = ac.aerlut.import_lut(lutid, lutdir, override=0)
            else:
                ## sensor specific lut
                lut_data_dict, lut_meta = ac.aerlut.import_lut_sensor(sensor, None, lutid, override=0, lutdir=lutdir)

                #bands = list(lut_data_dict.keys())
                # get bands from rsr_file as different systems may not keep dict keys in the same order
                rsr_file = ac.config['data_dir']+'/RSR/'+sensor+'.txt'
                rsr, rsr_bands = ac.shared.rsr_read(file=rsr_file)

            ## set up lut dimensions
            if ip == 0:
                # lut data is par, wave, azi, thv, ths, wind, tau
                lut_dim = [[i for i in pressures]] + [[i for i in range(len(lut_meta['par']))]]
                if sensor is None:
                    lutd = []
                    lut_dim += [lut_meta['wave']]

                lut_dim += [lut_meta[k] for k in ['azi', 'thv', 'ths', 'tau']]
                ipd = {par:ip for ip,par in enumerate(lut_meta['par'])}

                lut_dict[lut] = {'meta':lut_meta, 'dim':lut_dim, 'ipd':ipd}

                if sensor is not None:
                    lut_dict[lut]['lut'] = {band:[] for band in rsr_bands}

            if sensor is None:
                lutd.append(lut_data)
            else:
                for band in rsr_bands:
                    lut_dict[lut]['lut'][band].append(lut_data_dict[band])

        ## generic LUT
        if sensor is None:
            lut_dict[lut]['lut'] = np.stack(lutd)
            ipd = lut_dict[lut]['ipd']

            if add_rsky:
                tlut = lut_dict[lut]['lut']
                if add_rsky_winds:
                    if True:
                        rskyd = ac.aerlut.import_rsky_luts(models=[int(lut[-1])], lutbase=rsky_base, winds=rsky_winds)
                        rlut = rskyd[int(lut[-1])]['lut']
                        rskyd = None
                    else:
                        rlut = None
                        for iw, wind in enumerate(rsky_winds):
                            ret = ac.aerlut.import_rsky_lut(int(lut[-1]), lutbase=rsky_base.format(wind))
                            ret = {'lut':ret[0], 'meta':ret[1], 'dims':ret[2], 'rgi':ret[3]}
                            tmp = ret['lut']
                            rlut = tmp[:, :, :, :, None, :] if rlut is None else np.insert(rlut,(iw),tmp,axis=4)
                            ret = None
                            tmp = None
                    ## repeat 6sv lut for winds
                    tlut = np.repeat(tlut, len(rsky_winds), axis=-2)
                else:
                    ## import rsky lut
                    rlut, rmeta, rdim, rrgi = ac.aerlut.import_rsky_lut(int(lut[-1]), lutbase=rsky_lut)
                    ## add empty axis (wind placeholder)
                    rlut = rlut[:,:,:,:,None, :]

                ## repeat for pressures
                rlut = np.repeat(rlut[np.newaxis,:], len(pressures), axis=0)

                ## add to the LUT
                ## model rsky at surface
                ## rsky_s idx 22
                tlut = np.insert(tlut, (22), rlut, axis=1)

                ## model rsky at toa
                ## rsky_t idx 23
                ## (utott * dtott * rsky) / (1. - rsky * astot)
                tmp = (tlut[:, ipd['utott'],:,:,:,:,:]*\
                       tlut[:, ipd['dtott'],:,:,:,:,:]*
                       tlut[:, 22,:,:,:,:,:]) /\
                       (1.-tlut[:, 22,:,:,:,:,:] *\
                       tlut[:, ipd['astot'],:,:,:,:,:])
                tlut = np.insert(tlut, (23), tmp, axis=1)

                ## add romix+rsky
                ## idx 24
                tmp = tlut[:, ipd['romix'],:,:,:,:,:] + tlut[:, 23,:,:,:,:,:]
                tlut = np.insert(tlut, (24), tmp, axis=1)

                ## replace lut and add these parameters
                lut_dict[lut]['lut'] = tlut
                lut_dict[lut]['meta']['par'] += ['rsky_s', 'rsky_t', 'romix+rsky_t']
                lut_dict[lut]['ipd'] = {p:i for i,p in enumerate(lut_dict[lut]['meta']['par'])}
                lut_dict[lut]['dim'][1]+= [22, 23, 24]
                ## end add rsky

                if add_rsky_winds:
                    ## create new dim with winds
                    dim = [np.asarray(pressures),
                           np.asarray(lut_dict[lut]['dim'][1]),
                           lut_dict[lut]['meta']['wave'],
                           lut_dict[lut]['meta']['azi'],
                           lut_dict[lut]['meta']['thv'],
                           lut_dict[lut]['meta']['ths'],
                           np.asarray(rsky_winds),
                           lut_dict[lut]['meta']['tau']]
                    lut_dict[lut]['dim'] = dim
                ## end add rsky

            ## add product of transmittances
            if add_dutott:
                lut_dict[lut]['meta']['par']+=['dutott']
                ax = int(lut_dict[lut]['dim'][1][-1])+1
                lut_dict[lut]['dim'][1]=np.append(lut_dict[lut]['dim'][1],ax)
                lut_dict[lut]['ipd'] = {p:i for i,p in enumerate(lut_dict[lut]['meta']['par'])}
                iu = {p:i for i,p in enumerate(lut_dict[lut]['meta']['par'])}['utott']
                id = {p:i for i,p in enumerate(lut_dict[lut]['meta']['par'])}['dtott']
                tmp = lut_dict[lut]['lut'][:,iu,:,:,:,:,:,:]*lut_dict[lut]['lut'][:,id,:,:,:,:,:,:]
                lut_dict[lut]['lut'] = np.insert(lut_dict[lut]['lut'], (ax), tmp, axis=1)

            ## set up LUT interpolator
            if add_rsky & add_rsky_winds:
                lut_dict[lut]['rgi'] = scipy.interpolate.RegularGridInterpolator(lut_dict[lut]['dim'],
                                                                             lut_dict[lut]['lut'][:,:,:,:,:,:,:,:],
                                                                             bounds_error=False, fill_value=np.nan)
            else:
                lut_dict[lut]['rgi'] = scipy.interpolate.RegularGridInterpolator(lut_dict[lut]['dim'],
                                                                             lut_dict[lut]['lut'][:,:,:,:,:,:,0,:],
                                                                             bounds_error=False, fill_value=np.nan)
        else:
            ## make arrays
            for band in rsr_bands:
                lut_dict[lut]['lut'][band] = np.stack(lut_dict[lut]['lut'][band])

            ## add rsky if requested
            if add_rsky:
                if add_rsky_winds:
                    rskyd = ac.aerlut.import_rsky_luts(models=[int(lut[-1])], lutbase=rsky_base, winds=rsky_winds, sensor=sensor)
                    rlut = rskyd[int(lut[-1])]['lut']
                    rskyd = None
                else:
                    rlut, rmeta, rdim, rrgi = ac.aerlut.import_rsky_lut(int(lut[-1]), sensor=sensor, lutbase=rsky_lut)

                ## current pars
                ipd = {p:i for i,p in enumerate(lut_dict[lut]['meta']['par'])}

                ## run through bands
                for band in rlut:
                    #if rsky_lut[0:23] == 'ACOLITE-RSKY-202101-75W':
                    #    tmp = rlut[band] * 1.0
                    tmp = rlut[band] * 1.0
                    tlut = lut_dict[lut]['lut'][band]
                    if add_rsky_winds:
                        tmp = rlut[band] * 1.0
                        tlut = np.repeat(tlut, len(rsky_winds), axis=-2)
                    else:
                        ## add empty axis (wind placeholder)
                        tmp = tmp[:,:,:,None, :]

                    ## add to the LUT
                    ## model rsky at surface
                    ## rsky_s idx 22
                    tlut = np.insert(tlut, (22), tmp, axis=1)

                    ## model rsky at toa
                    ## rsky_t idx 23
                    ## (utott * dtott * rsky) / (1. - rsky * astot)
                    tmp = (tlut[:, ipd['utott'],:,:,:,:,:]*\
                           tlut[:, ipd['dtott'],:,:,:,:,:]*
                           tlut[:, 22,:,:,:,:,:]) /\
                          (1.-tlut[:, 22,:,:,:,:,:] *\
                           tlut[:, ipd['astot'],:,:,:,:,:])
                    tlut = np.insert(tlut, (23), tmp, axis=1)

                    ## add romix+rsky
                    ## idx 24
                    tmp = tlut[:, ipd['romix'],:,:,:,:,:] + tlut[:, 23,:,:,:,:,:]
                    tlut = np.insert(tlut, (24), tmp, axis=1)

                    ## replace in dict
                    lut_dict[lut]['lut'][band] = tlut

                ## add new pars
                lut_dict[lut]['meta']['par'] += ['rsky_s', 'rsky_t', 'romix+rsky_t']
                lut_dict[lut]['dim'][1]+= [22, 23, 24]
                if add_rsky_winds:
                    ## create new dim with winds
                    dim = [np.asarray(pressures),
                           np.asarray(lut_dict[lut]['dim'][1]),
                           lut_dict[lut]['meta']['azi'],
                           lut_dict[lut]['meta']['thv'],
                           lut_dict[lut]['meta']['ths'],
                           np.asarray(rsky_winds),
                           lut_dict[lut]['meta']['tau']]
                    lut_dict[lut]['dim'] = dim
                ## end add rsky

            ## add product of transmittances
            if add_dutott:
                lut_dict[lut]['meta']['par']+=['dutott']
                ax = int(lut_dict[lut]['dim'][1][-1])+1
                lut_dict[lut]['dim'][1]=np.append(lut_dict[lut]['dim'][1],ax)
                iu = {p:i for i,p in enumerate(lut_dict[lut]['meta']['par'])}['utott']
                id = {p:i for i,p in enumerate(lut_dict[lut]['meta']['par'])}['dtott']
                for band in rsr_bands:
                    tmp = lut_dict[lut]['lut'][band][:,iu,:,:,:,:,:]*lut_dict[lut]['lut'][band][:,id,:,:,:,:,:]
                    lut_dict[lut]['lut'][band] = np.insert(lut_dict[lut]['lut'][band], (ax), tmp, axis=1)

            lut_dict[lut]['rgi'] = {}
            lut_dict[lut]['ipd'] = {p:i for i,p in enumerate(lut_dict[lut]['meta']['par'])}
            for band in rsr_bands:
                ## set up LUT interpolator per band
                if add_rsky & add_rsky_winds:
                    lut_dict[lut]['rgi'][band] = scipy.interpolate.RegularGridInterpolator(lut_dict[lut]['dim'],
                                                                                       lut_dict[lut]['lut'][band][:,:,:,:,:,:,:],
                                                                                       bounds_error=False, fill_value=np.nan)
                else:
                    lut_dict[lut]['rgi'][band] = scipy.interpolate.RegularGridInterpolator(lut_dict[lut]['dim'],
                                                                                       lut_dict[lut]['lut'][band][:,:,:,:,:,0,:],
                                                                                       bounds_error=False, fill_value=np.nan)

    return(lut_dict)
