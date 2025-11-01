## read all luts and set up rgi
## QV 2020-01-15
## Last modifications: 2020-07-14 (QV) added sensor option
##                     2020-07-14 (QV) changed sensor option to a dict per band
##                     2020-07-22 (QV) added the rsky lut to the lut and rgi - need to update rsky azimuths
##                     2020-07-25 (QV) added updates to rsky lut
##                     2021-01-19 (QV) update to new LUTs
##                     2021-02-01 (QV) added wind speed for rsky lut, removed temp fixes for old luts
##                     2021-02-03 (QV) added down*up total transmittances
##                     2021-03-01 (QV) added new rsky luts with integrated wind speed
##                     2021-06-08 (QV) added lut par subsetting
##                     2021-10-24 (QV) added get_remote as keyword
##                     2021-11-09 (QV) added reduce dimensions
##                     2022-03-03 (QV) increased default reduce dimensions AOT range
##                     2024-04-25 (QV) check par/add_rsky combination, check wnd dimensions
##                                     separate function for merging luts - which is not faster (maybe depending on disk speed?)

def import_luts(pressures = [500, 750, 1013, 1100],
                base_luts = ['ACOLITE-LUT-202110-MOD1', 'ACOLITE-LUT-202110-MOD2'],
                rsky_lut = 'ACOLITE-RSKY-202102-82W',
                lut_par = ['utott', 'dtott', 'astot', 'ttot', 'romix'],
                reduce_dimensions = False, return_lut_array = False,
                par = 'romix',
                vza_range = [0, 16],  aot_range = [0, 1.5],
                use_merged_lut = False, store_merged_lut = False,
                get_remote = True, sensor = None, add_rsky = False, add_dutott = True):
    import scipy.interpolate
    import numpy as np
    import acolite as ac

    ## indices for reducing LUT size
    vza_sub = [0, -1]
    aot_sub = [0, -1]

    if par not in ['romix+rsky_t', 'romix+ffss_toa', 'romix+rsurf']: add_rsky = False

    if add_rsky:
        klist = ['utott', 'dtott', 'astot', 'romix']
        if (par == 'romix+rsky_t') & (use_merged_lut): klist += ['romix+rsky_t']
        if (par == 'romix+ffss_toa') & (use_merged_lut): klist += ['romix+ffss_toa']
        if par == 'romix+rsurf': klist += ['rsurf']

        for k in klist:
            if k not in lut_par: lut_par.append(k)
    if add_dutott:
        for k in ['utott', 'dtott']:
            if k not in lut_par: lut_par.append(k)

    lut_dict = {}
    ## run through luts
    for lut in base_luts:
        ## run through pressures
        for ip, pr in enumerate(pressures):
            lutid = '{}-{}mb'.format(lut, '{}'.format(pr).zfill(4))
            lutdir = '{}/{}'.format(ac.config['lut_dir'], '-'.join(lutid.split('-')[0:3]))

            if sensor is None:
                ## indices for reducing LUT size
                vza_idx, aot_idx = 4, 7
                if not add_rsky: aot_idx = 6
                if (add_rsky) & (par == 'romix+rsky_t') & (use_merged_lut): ## load merged lut
                    lut_data, lut_meta = ac.aerlut.merged_lut(lut, rsky_lut, pr, sensor = sensor,
                                                              store = store_merged_lut, lut_par = lut_par, get_remote = get_remote)
                else:
                    lut_data, lut_meta = ac.aerlut.import_lut(lutid, lutdir, sensor = sensor, lut_par = lut_par, get_remote = get_remote)
            else:
                ## indices for reducing LUT size
                vza_idx, aot_idx = 3, 6
                if not add_rsky: aot_idx = 5

                if (add_rsky) & (par == 'romix+rsky_t') & (use_merged_lut):  ## load merged lut
                    lut_data_dict, lut_meta = ac.aerlut.merged_lut(lut, rsky_lut, pr, sensor = sensor,
                                                                   store = store_merged_lut, lut_par = lut_par, get_remote = get_remote)
                else:
                    lut_data_dict, lut_meta = ac.aerlut.import_lut(lutid, lutdir, sensor = sensor, lut_par = lut_par, get_remote = get_remote)

                if 'bands' not in lut_meta:
                    # get bands from rsr_file as different systems may not keep dict keys in the same order
                    rsr_file = ac.config['data_dir']+'/RSR/'+sensor+'.txt'
                    rsr, rsr_bands = ac.shared.rsr_read(file=rsr_file)
                else:
                    rsr_bands = lut_meta['bands']

            ## set up lut dimensions
            if ip == 0:
                # lut data is par, wave, azi, thv, ths, wind, tau
                lut_dim = [[i for i in pressures]] + [[i for i in range(len(lut_meta['par']))]]
                if sensor is None:
                    lutd = []
                    lut_dim += [lut_meta['wave']]
                if (add_rsky) & ((par == 'romix+rsky_t') | (par == 'romix+ffss_toa')) & (use_merged_lut): ## add wind dimension
                    lut_dim += [lut_meta[k] for k in ['azi', 'thv', 'ths', 'wnd', 'tau']]
                else:
                    lut_dim += [lut_meta[k] for k in ['azi', 'thv', 'ths', 'tau']]
                ipd = {par:ip for ip,par in enumerate(lut_meta['par'])}

                lut_dict[lut] = {'meta':lut_meta, 'dim':lut_dim, 'ipd':ipd}

                ## find indices to reduce dimensions
                if reduce_dimensions:
                    for vi, v in enumerate(lut_meta['thv']):
                        if v <= vza_range[0]: vza_sub[0] = vi
                        if (v >= vza_range[1]) & (vza_sub[1] == -1): vza_sub[1] = vi+1
                    for vi, v in enumerate(lut_meta['tau']):
                        if v <= aot_range[0]: aot_sub[0] = vi
                        if (v >= aot_range[1]) & (aot_sub[1] == -1): aot_sub[1] = vi+1

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

            if (add_rsky) & ((par == 'romix+rsky_t') | (par == 'romix+ffss_toa')) & (not use_merged_lut):
                tlut = lut_dict[lut]['lut']

                if (par == 'romix+rsky_t'):
                    rskyd = ac.aerlut.import_rsky_luts(models=[int(lut[-1])], lutbase=rsky_lut, get_remote = get_remote)
                    rlut = rskyd[int(lut[-1])]['lut']
                    rsky_winds  = rskyd[int(lut[-1])]['meta']['wind']
                    rskyd = None
                    add_par = ['rsky_s', 'rsky_t', 'romix+rsky_t']
                if (par == 'romix+ffss_toa'):
                    meta_rsky, rlut = ac.ac.skydome.import_skydome_lut_reformat(sensor = None, model = int(lut[-1])-1)
                    rsky_winds  = meta_rsky['wind']
                    add_par = ['ffss_boa', 'ffss_toa', 'romix+ffss_toa']

                ## repeat 6sv lut for winds
                tlut = np.repeat(tlut, len(rsky_winds), axis=-2)

                ## repeat for pressures
                rlut = np.repeat(rlut[np.newaxis,:], len(pressures), axis=0)

                ## add to the LUT
                ## model rsky at surface
                ax = len(ipd)
                tlut = np.insert(tlut, (ax), rlut, axis=1)

                ## model rsky at toa
                ## (utott * dtott * rsky) / (1. - rsky * astot)
                tmp = (tlut[:, ipd['utott'],:,:,:,:,:]*\
                       tlut[:, ipd['dtott'],:,:,:,:,:]*
                       tlut[:, ax,:,:,:,:,:]) /\
                       (1.-tlut[:, ax,:,:,:,:,:] *\
                       tlut[:, ipd['astot'],:,:,:,:,:])
                tlut = np.insert(tlut, (ax+1), tmp, axis=1)

                ## add romix+rsky
                tmp = tlut[:, ipd['romix'],:,:,:,:,:] + tlut[:, ax+1,:,:,:,:,:]
                tlut = np.insert(tlut, (ax+2), tmp, axis=1)

                ## replace lut and add these parameters
                lut_dict[lut]['lut'] = tlut
                lut_dict[lut]['meta']['par'] += add_par
                lut_dict[lut]['ipd'] = {p:i for i,p in enumerate(lut_dict[lut]['meta']['par'])}
                lut_dict[lut]['dim'][1]+= [ax, ax+1, ax+2]
                ## end add rsky

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
                ## end add rsky romix+rsky_t

            if (add_rsky) & (par == 'romix+rsurf'):
                ax = len(ipd)
                tmp = lut_dict[lut]['lut'][:, ipd['romix'],:,:,:,:,:] + lut_dict[lut]['lut'][:, ipd['rsurf'],:,:,:,:,:]
                lut_dict[lut]['lut'] = np.insert(lut_dict[lut]['lut'], (ax), tmp, axis=1)
                tmp = None
                lut_dict[lut]['meta']['par'] += ['romix+rsurf']
                lut_dict[lut]['ipd'] = {p:i for i,p in enumerate(lut_dict[lut]['meta']['par'])}
                lut_dict[lut]['dim'][1]+= [ax]
                ## create new dim with winds
                dim = [np.asarray(pressures),
                       np.asarray(lut_dict[lut]['dim'][1]),
                       lut_dict[lut]['meta']['wave'],
                       lut_dict[lut]['meta']['azi'],
                       lut_dict[lut]['meta']['thv'],
                       lut_dict[lut]['meta']['ths'],
                       np.atleast_1d(lut_dict[lut]['meta']['wnd']),
                       lut_dict[lut]['meta']['tau']]
                lut_dict[lut]['dim'] = dim

            ## reduce dimensions / memory use
            if reduce_dimensions:
                lut_dict[lut]['dim'][vza_idx] = lut_dict[lut]['dim'][vza_idx][vza_sub[0]:vza_sub[1]]
                lut_dict[lut]['dim'][aot_idx] = lut_dict[lut]['dim'][aot_idx][aot_sub[0]:aot_sub[1]]
                lut_dict[lut]['lut'] = lut_dict[lut]['lut'][:,:,:,:,vza_sub[0]:vza_sub[1],:,:,aot_sub[0]:aot_sub[1]]

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
            if add_rsky:
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
            if (add_rsky) & ((par == 'romix+rsky_t') | (par == 'romix+ffss_toa')) & (not use_merged_lut):
                expand_dims = False
                if (par == 'romix+rsky_t'):
                    rskyd = ac.aerlut.import_rsky_luts(models = [int(lut[-1])], lutbase = rsky_lut, sensor = sensor, get_remote = get_remote)
                    rlut = rskyd[int(lut[-1])]['lut']
                    if 'wind' in rskyd[int(lut[-1])]['meta']:
                        rsky_winds  = rskyd[int(lut[-1])]['meta']['wind']
                    else:
                        rsky_winds = np.atleast_1d([0,20])
                        expand_dims = True
                    add_par = ['rsky_s', 'rsky_t', 'romix+rsky_t']
                if (par == 'romix+ffss_toa'):
                    meta_rsky, rlut = ac.ac.skydome.import_skydome_lut_reformat(sensor = sensor, model = int(lut[-1])-1)
                    rsky_winds  = meta_rsky['wind']
                    add_par = ['ffss_boa', 'ffss_toa', 'romix+ffss_toa']

                ## current pars
                ipd = {p:i for i,p in enumerate(lut_dict[lut]['meta']['par'])}

                ## run through bands
                for band in rlut:
                    tmp = rlut[band] * 1.0
                    tlut = lut_dict[lut]['lut'][band]
                    tmp = rlut[band] * 1.0

                    ## add wind axis
                    if expand_dims:
                        tmp = np.expand_dims(tmp, -2)
                        tmp = np.repeat(tmp, len(rsky_winds), axis=-2)
                    tlut = np.repeat(tlut, len(rsky_winds), axis=-2)

                    ## add to the LUT
                    ## model rsky at surface
                    ax = len(ipd)
                    tlut = np.insert(tlut, (ax), tmp, axis=1)

                    ## model rsky at toa
                    ## (utott * dtott * rsky) / (1. - rsky * astot)
                    tmp = (tlut[:, ipd['utott'],:,:,:,:,:]*\
                           tlut[:, ipd['dtott'],:,:,:,:,:]*
                           tlut[:, ax,:,:,:,:,:]) /\
                          (1.-tlut[:, ax,:,:,:,:,:] *\
                           tlut[:, ipd['astot'],:,:,:,:,:])
                    tlut = np.insert(tlut, (ax+1), tmp, axis=1)

                    ## add romix+rsky
                    tmp = tlut[:, ipd['romix'],:,:,:,:,:] + tlut[:, ax+1,:,:,:,:,:]
                    tlut = np.insert(tlut, (ax+2), tmp, axis=1)

                    ## replace in dict
                    lut_dict[lut]['lut'][band] = tlut

                ## add new pars
                lut_dict[lut]['meta']['par'] += add_par
                lut_dict[lut]['dim'][1]+= [ax, ax+1, ax+2]

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

            if (add_rsky) & (par == 'romix+rsurf'):
                ## current pars
                ipd = {p:i for i,p in enumerate(lut_dict[lut]['meta']['par'])}
                ax = len(ipd)

                ## run through bands
                for band in lut_dict[lut]['lut']:
                    tmp = lut_dict[lut]['lut'][band][:, ipd['romix'],:,:,:,:] + lut_dict[lut]['lut'][band][:, ipd['rsurf'],:,:,:,:]
                    lut_dict[lut]['lut'][band] = np.insert(lut_dict[lut]['lut'][band], (ax), tmp, axis=1)
                    tmp = None
                lut_dict[lut]['meta']['par'] += ['romix+rsurf']
                lut_dict[lut]['ipd'] = {p:i for i,p in enumerate(lut_dict[lut]['meta']['par'])}
                lut_dict[lut]['dim'][1]+= [ax]
                ## create new dim with winds
                dim = [np.asarray(pressures),
                       np.asarray(lut_dict[lut]['dim'][1]),
                       lut_dict[lut]['meta']['azi'],
                       lut_dict[lut]['meta']['thv'],
                       lut_dict[lut]['meta']['ths'],
                       np.atleast_1d(lut_dict[lut]['meta']['wnd']),
                       lut_dict[lut]['meta']['tau']]
                lut_dict[lut]['dim'] = dim

            ## reduce dimensions / memory use
            if reduce_dimensions:
                lut_dict[lut]['dim'][vza_idx] = lut_dict[lut]['dim'][vza_idx][vza_sub[0]:vza_sub[1]]
                lut_dict[lut]['dim'][aot_idx] = lut_dict[lut]['dim'][aot_idx][aot_sub[0]:aot_sub[1]]
                for band in rsr_bands:
                    lut_dict[lut]['lut'][band] = lut_dict[lut]['lut'][band][:,:,:,vza_sub[0]:vza_sub[1],:,:,aot_sub[0]:aot_sub[1]]

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
                if add_rsky:
                    lut_dict[lut]['rgi'][band] = scipy.interpolate.RegularGridInterpolator(lut_dict[lut]['dim'],
                                                                                       lut_dict[lut]['lut'][band][:,:,:,:,:,:,:],
                                                                                       bounds_error=False, fill_value=np.nan)
                else:
                    lut_dict[lut]['rgi'][band] = scipy.interpolate.RegularGridInterpolator(lut_dict[lut]['dim'],
                                                                                       lut_dict[lut]['lut'][band][:,:,:,:,:,0,:],
                                                                                       bounds_error=False, fill_value=np.nan)
    ## remove LUT array to reduce memory use
    if not return_lut_array:
        for lut in lut_dict:
            del lut_dict[lut]['lut']

    return(lut_dict)
