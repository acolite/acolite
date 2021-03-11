## def acolite_gem
## runs ACOLITE/DSF on generic extracted miniscene
## written by Quinten Vanhellemont, RBINS
## 2021-03-01
## modifications:


def acolite_gem(gem,
                output = None,
                sub = None,
                settings = None,

                output_file = True,
                target_file = None,

                return_gem = False,

                verbosity=0):

    import os, time
    import numpy as np
    import scipy.ndimage
    import acolite as ac

    ## for gaussian smoothing of aot
    from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
    from astropy.convolution import convolve

    ## read gem file if NetCDF
    if type(gem) is str:
        gemf = '{}'.format(gem)
        gem = ac.gem.read(gemf, sub=sub)
    gemf = gem['gatts']['gemfile']

    ## combine default and user defined settings
    setu = ac.acolite.settings.parse(gem['gatts']['sensor'], settings=settings)

    ## check blackfill
    if setu['blackfill_skip']:
        rhot_ds = [ds for ds in gem['datasets'] if 'rhot_' in ds]
        rhot_wv = [int(ds.split('_')[1]) for ds in rhot_ds]
        bi, bw = ac.shared.closest_idx(rhot_wv, setu['blackfill_wave'])
        #band_data = 1.0*gem['data'][rhot_ds[bi]]
        npx = gem['data'][rhot_ds[bi]].shape[0] * gem['data'][rhot_ds[bi]].shape[1]
        #nbf = npx - len(np.where(np.isfinite(gem['data'][rhot_ds[bi]])*(gem['data'][rhot_ds[bi]]>0))[0])
        nbf = npx - len(np.where(np.isfinite(gem['data'][rhot_ds[bi]]))[0])
        #band_data = None
        if (nbf/npx) >= float(setu['blackfill_max']):
            if verbosity>1: print('Skipping scene as crop is {:.0f}% blackfill'.format(100*nbf/npx))
            return()

    if verbosity > 0: print('Running acolite for {}'.format(gemf))

    output_name = gem['gatts']['output_name'] if 'output_name' in gem['gatts'] else os.path.basename(gemf).replace('.nc', '')

    ## get dimensions and number of elements
    gem['gatts']['data_dimensions'] = gem['data']['lon'].shape
    gem['gatts']['data_elements'] = gem['gatts']['data_dimensions'][0]*gem['gatts']['data_dimensions'][1]

    ## read rsrd and get band wavelengths
    rsrd = ac.shared.rsr_dict(gem['gatts']['sensor'])
    if gem['gatts']['sensor'] in rsrd:
        rsrd = rsrd[gem['gatts']['sensor']]
    else:
        rsrd = None
        stop

    ## set defaults
    gem['gatts']['uoz'] = setu['uoz_default']
    gem['gatts']['uwv'] = setu['uwv_default']
    gem['gatts']['wind'] = setu['wind']
    gem['gatts']['pressure'] = setu['pressure']

    ## read ancillary data
    if setu['ancillary_data']:
        clon = np.nanmedian(gem['data']['lon'])
        clat = np.nanmedian(gem['data']['lat'])
        anc = ac.ac.ancillary.get(gem['gatts']['isodate'], clon, clat)

        ## overwrite the defaults
        if ('ozone' in anc): gem['gatts']['uoz'] = anc['ozone']['interp']/1000. ## convert from MET data
        if ('p_water' in anc): gem['gatts']['uwv'] = anc['p_water']['interp']/10. ## convert from MET data
        if ('z_wind' in anc) & ('m_wind' in anc):
            gem['gatts']['wind'] = ((anc['z_wind']['interp']**2) + (anc['m_wind']['interp']**2))**0.5
        if ('press' in anc): gem['gatts']['pressure'] = anc['press']['interp']

    ## get mean average geometry
    geom_ds = ['sza', 'vza', 'raa', 'pressure', 'wind']
    geom_mean = {k: np.nanmean(gem['data'][k]) if k in gem['data'] else gem['gatts'][k] for k in geom_ds}

    ## set wind to wind range
    gem['gatts']['wind'] = max(0.1, gem['gatts']['wind'])
    gem['gatts']['wind'] = min(20, gem['gatts']['wind'])

    ## get gas transmittance
    tg_dict = ac.ac.gas_transmittance(geom_mean['sza'], geom_mean['vza'],
                                      uoz=gem['gatts']['uoz'], uwv=gem['gatts']['uwv'],
                                      sensor=gem['gatts']['sensor'])

    ## make bands dataset
    gem['bands'] = {}
    for bi, b in enumerate(rsrd['rsr_bands']):
        if b not in gem['bands']:
            gem['bands'][b] = {k:rsrd[k][b] for k in ['wave_mu', 'wave_nm', 'wave_name'] if b in rsrd[k]}
            gem['bands'][b]['rhot_ds'] = 'rhot_{}'.format(gem['bands'][b]['wave_name'])
            gem['bands'][b]['rhos_ds'] = 'rhos_{}'.format(gem['bands'][b]['wave_name'])
            for k in tg_dict:
                if k not in ['wave']: gem['bands'][b][k] = tg_dict[k][b]
    ## end bands dataset

    ## determine use of reverse lut rhot->aot
    use_revlut = False
    ## if path reflectance is tiled or resolved, use reverse lut
    #if setu['dsf_path_reflectance'] != 'fixed': use_revlut = True
    ## no need to use reverse lut if fixed geometry is used
    #if setu['resolved_geometry']:
    ## we want to use the reverse lut to derive aot if the geometry data is resolved
    for ds in geom_ds:
        if ds not in gem['data']:
            gem['data'][ds] = geom_mean[ds]
        if len(np.atleast_1d(gem['data'][ds]))>1:
            use_revlut=True ## if any dataset more than 1 dimension use revlut
        gem['data']['{}_mean'.format(ds)] = np.asarray(np.nanmean(gem['data'][ds])) ## also store tile mean
        gem['data']['{}_mean'.format(ds)].shape+=(1,1) ## make 1,1 dimensions
    if not setu['resolved_geometry']: use_revlut = False

    ## for ease of subsetting later, repeat single element datasets to the tile shape
    if use_revlut:
        for ds in geom_ds:
            if len(np.atleast_1d(gem['data'][ds]))!=1: continue
            gem['data'][ds] = np.repeat(gem['data'][ds], gem['gatts']['data_elements']).reshape(gem['gatts']['data_dimensions'])
    else:
        ## if use revlut is False at this point, we don't have resolved geometry
        setu['resolved_geometry'] = False
    #print(use_revlut, setu['dsf_path_reflectance'], setu['resolved_geometry'])
    ## end determine revlut

    ## for tiled processing track tile positions and average geometry
    tiles = []
    if (setu['dsf_path_reflectance'] == 'tiled') & (setu['dsf_tile_dimensions'] is not None):
        ni = np.ceil(gem['gatts']['data_dimensions'][0]/setu['dsf_tile_dimensions'][0]).astype(int)
        nj = np.ceil(gem['gatts']['data_dimensions'][1]/setu['dsf_tile_dimensions'][1]).astype(int)
        if (ni <= 1) | (nj <= 1):
            if verbosity > 1: print('Scene too small for tiling ({}x{} tiles), using fixed processing'.format(ni,nj))
            setu['dsf_path_reflectance'] = 'fixed'
        else:
            ntiles = ni*nj
            if verbosity > 1: print('Processing with {} tiles ({}x{})'.format(ntiles, ni, nj))

            ## compute tile dimensions
            for ti in range(ni):
                for tj in range(nj):
                    subti = [setu['dsf_tile_dimensions'][0]*ti, setu['dsf_tile_dimensions'][0]*(ti+1)]
                    subti[1] = np.min((subti[1], gem['gatts']['data_dimensions'][0]))
                    subtj = [setu['dsf_tile_dimensions'][1]*tj, setu['dsf_tile_dimensions'][1]*(tj+1)]
                    subtj[1] = np.min((subtj[1], gem['gatts']['data_dimensions'][1]))
                    tiles.append((ti, tj, subti, subtj))

            ## create tile geometry datasets
            for ds in geom_ds:
                if len(np.atleast_1d(gem['data'][ds]))>1: ## if not fixed geometry
                    gem['data']['{}_tiled'.format(ds)] = np.zeros((ni,nj))+np.nan
                    for t in range(ntiles):
                        ti, tj, subti, subtj = tiles[t]
                        gem['data']['{}_tiled'.format(ds)][ti, tj] = \
                            np.nanmean(gem['data'][ds][subti[0]:subti[1],subtj[0]:subtj[1]])
                else: ## if fixed geometry
                    gem['data']['{}_tiled'.format(ds)] = gem['data'][ds]
    ## end tiling

    ## read LUTs
    if setu['sky_correction']:
        par = 'romix+rsky_t'
    else:
        par = 'romix'

    ## load reverse lut romix -> aot
    if use_revlut: revl = ac.aerlut.reverse_lut(gem['gatts']['sensor'], par=par)
    ## load aot -> atmospheric parameters lut
    lutdw = ac.aerlut.import_luts(add_rsky=True, sensor=gem['gatts']['sensor'])
    luts = list(lutdw.keys())

    ## run through bands to get aot
    aot_dict = {}
    dsf_rhod = {}
    for bi, b in enumerate(gem['bands']):
        if (b in setu['dsf_exclude_bands']): continue
        if ('rhot_ds' not in gem['bands'][b]) or ('tt_gas' not in gem['bands'][b]): continue
        if gem['bands'][b]['rhot_ds'] not in gem['data']: continue

        ## skip band for aot computation
        if gem['bands'][b]['tt_gas'] < setu['min_tgas_aot']: continue

        ## skip bands according to configuration
        if (gem['bands'][b]['wave_nm'] < setu['dsf_wave_range'][0]): continue
        if (gem['bands'][b]['wave_nm'] > setu['dsf_wave_range'][1]): continue
        if (b in setu['dsf_exclude_bands']): continue

        if verbosity > 1: print(b, gem['bands'][b]['rhot_ds'])

        band_data = gem['data'][gem['bands'][b]['rhot_ds']]*1.0
        valid = np.isfinite(band_data)*(band_data>0)
        mask = valid is False

        ## apply TOA filter
        if setu['dsf_filter_toa']:
            band_data[mask] = np.nanmedian(band_data) ## fill mask with median
            band_data = scipy.ndimage.median_filter(band_data, size=setu['dsf_filter_box'])
            band_data = scipy.ndimage.percentile_filter(band_data, setu['dsf_filter_percentile'], size=setu['dsf_filter_box'])
            band_data[mask] = np.nan
        band_sub = np.where(valid)

        ## geometry key '' if using resolved, otherwise '_mean' or '_tiled'
        gk = ''

        ## fixed path reflectance
        if setu['dsf_path_reflectance'] == 'fixed':
            if setu['dsf_spectrum_option'] == 'darkest':
                band_data = np.array((np.nanpercentile(band_data[band_sub], 0)))
            if setu['dsf_spectrum_option'] == 'percentile':
                band_data = np.array((np.nanpercentile(band_data[band_sub], setu['dsf_percentile'])))
            if setu['dsf_spectrum_option'] == 'intercept':
                band_data = ac.shared.intercept(band_data[band_sub], setu['dsf_intercept_pixels'])
            band_data.shape+=(1,1) ## make 1,1 dimensions
            if verbosity > 2: print(b, setu['dsf_spectrum_option'], '{:.3f}'.format(band_data[0,0]))
            gk='_mean'

        ## tiled path reflectance
        elif setu['dsf_path_reflectance'] == 'tiled':
            gk = '_tiled'

            ## tile this band data
            tile_data = np.zeros((tiles[-1][0]+1, tiles[-1][1]+1)) + np.nan
            for t in range(len(tiles)):
                ti, tj, subti, subtj = tiles[t]
                tsub = band_data[subtj[0]:subtj[1], subti[0]:subti[1]]
                tel = (subtj[1]-subtj[0]) * (subti[1]-subti[0])
                nsub = len(np.where(np.isfinite(tsub))[0])
                if nsub < tel * float(setu['dsf_min_tile_cover']): continue

                ## get per tile darkest
                if setu['dsf_spectrum_option'] == 'darkest':
                    tile_data[ti,tj] = np.array((np.nanpercentile(tsub, 0)))
                if setu['dsf_spectrum_option'] == 'percentile':
                    tile_data[ti,tj] = np.array((np.nanpercentile(tsub, setu['dsf_percentile'])))
                if setu['dsf_spectrum_option'] == 'intercept':
                    tile_data[ti,tj] = ac.shared.intercept(tsub, int(setu['dsf_intercept_pixels']))

            ## fill nan tiles with closest values
            ind = scipy.ndimage.distance_transform_edt(np.isnan(tile_data), return_distances=False, return_indices=True)
            band_data = tile_data[tuple(ind)]

        ## resolved per pixel dsf
        elif setu['dsf_path_reflectance'] == 'resolved':
            if not setu['resolved_geometry']: gk = '_mean'
        else:
            print('DSF option {} not configured'.format(setu['dsf_path_reflectance']))
            continue

        ## do gas correction
        band_sub = np.where(np.isfinite(band_data))
        band_data[band_sub] /= gem['bands'][b]['tt_gas']

        ## store rhod
        if setu['dsf_path_reflectance'] in ['fixed', 'tiled']:
            dsf_rhod[b] = band_data

        ## compute aot
        aot_band = {}
        for li, lut in enumerate(luts):
            aot_band[lut] = np.zeros(band_data.shape)+np.nan
            t0 = time.time()

            ## reverse lut interpolates rhot directly to aot
            if use_revlut:
                aot_band[lut][band_sub] = revl[lut]['rgi'][b]((gem['data']['pressure'+gk][band_sub],
                                                               gem['data']['raa'+gk][band_sub],
                                                               gem['data']['vza'+gk][band_sub],
                                                               gem['data']['sza'+gk][band_sub],
                                                               gem['data']['wind'+gk][band_sub],
                                                               band_data[band_sub]))
                # mask out of range aot
                aot_band[lut][aot_band[lut]<=revl[lut]['minaot']]=revl[lut]['minaot']#np.nan
                aot_band[lut][aot_band[lut]>=revl[lut]['maxaot']]=revl[lut]['maxaot']#np.nan

            ## standard lut interpolates rhot to results for different aot values
            else:
                ## get rho path for lut steps in aot
                tmp = lutdw[lut]['rgi'][b]((gem['data']['pressure'+gk],
                                            lutdw[lut]['ipd'][par],
                                            gem['data']['raa'+gk],
                                            gem['data']['vza'+gk],
                                            gem['data']['sza'+gk],
                                            gem['data']['wind'+gk], lutdw[lut]['meta']['tau']))
                tmp = tmp.flatten()
                #if setu['dsf_path_reflectance'] == 'fixed':
                #    print(gem['data']['pressure'+gk],gem['data']['raa'+gk],
                #                gem['data']['vza'+gk], gem['data']['sza'+gk],gem['data']['wind'+gk])

                ## interpolate rho path to observation
                aot_band[lut][band_sub] = np.interp(band_data[band_sub], tmp,
                                                   lutdw[lut]['meta']['tau'],
                                                   left=lutdw[lut]['meta']['tau'][0],
                                                   right=lutdw[lut]['meta']['tau'][-1])#left=np.nan, right=np.nan)

            tel = time.time()-t0

            if verbosity > 1: print('{}/B{} {} took {:.3f}s ({})'.format(gem['gatts']['sensor'], b, lut, tel, 'RevLUT' if use_revlut else 'StdLUT'))

        ###
        aot_dict[b] = aot_band

    ## get and sort keys
    aot_bands = list(aot_dict.keys())
    aot_bands.sort()

    ## get min aot per pixel
    aot_stack = {}
    for li, lut in enumerate(luts):
        ## stack aot for this lut
        for bi, b in enumerate(aot_bands):
            if b not in aot_dict: continue
            if lut not in aot_stack:
                aot_stack[lut] = {'all':  aot_dict[b][lut]*1.0}
            else:
                aot_stack[lut]['all'] = np.dstack((aot_stack[lut]['all'],
                                                   aot_dict[b][lut]))
        ## get minimum and mask of aot
        aot_stack[lut]['min'] = np.nanmin(aot_stack[lut]['all'], axis=2)
        aot_stack[lut]['mask'] = ~np.isfinite(aot_stack[lut]['min'])


        ## apply percentile filter
        if (setu['dsf_filter_aot']) & (setu['dsf_path_reflectance'] == 'resolved'):
            aot_stack[lut]['min'] = \
                scipy.ndimage.percentile_filter(aot_stack[lut]['min'],
                                                setu['dsf_filter_percentile'],
                                                size=setu['dsf_filter_box'])
        ## apply gaussian kernel smoothing
        if (setu['dsf_smooth_aot']) & (setu['dsf_path_reflectance'] == 'resolved'):
            aot_stack[lut]['min'] = \
                convolve(aot_stack[lut]['min'],
                         Gaussian2DKernel(x_stddev=setu['dsf_smooth_box'][0], y_stddev=setu['dsf_smooth_box'][1]),
                         boundary='extend')

        ## mask aot
        aot_stack[lut]['min'][aot_stack[lut]['mask']] = np.nan

        ## fill nan tiles with closest values
        #if (setu['dsf_path_reflectance'] == 'tiled'):
        #    ind = scipy.ndimage.distance_transform_edt(aot_stack[lut]['mask'],
        #                                               return_distances=False, return_indices=True)
        #    aot_stack[lut]['min'] = aot_stack[lut]['min'][tuple(ind)]

        ## store b1 and b2
        tmp = np.argsort(aot_stack[lut]['all'], axis=2)
        aot_stack[lut]['b1'] = tmp[:,:,0].astype(int)#.astype(float)
        #aot_stack[lut]['b1'][aot_stack[lut]['mask']] = np.nan
        aot_stack[lut]['b1'][aot_stack[lut]['mask']] = -1

        aot_stack[lut]['b2'] = tmp[:,:,1].astype(int)#.astype(float)
        #aot_stack[lut]['b2'][aot_stack[lut]['mask']] = np.nan
        aot_stack[lut]['b2'][aot_stack[lut]['mask']] = -1

        if setu['dsf_model_selection'] == 'min_dtau':
            ## array idices
            aid = np.indices(aot_stack[lut]['all'].shape[0:2])
            ## abs difference between first and second band tau
            aot_stack[lut]['dtau'] = np.abs(aot_stack[lut]['all'][aid[0,:],aid[1,:],tmp[:,:,0]]-\
                                            aot_stack[lut]['all'][aid[0,:],aid[1,:],tmp[:,:,1]])
            #return(aot_stack)
            #aot_stack[lut]['dtau'] =
        tmp = None

    ## select model based on min rmsd for 2 bands
    if verbosity > 1: print('Choosing best fitting model: {}'.format(setu['dsf_model_selection']))

    ## run through model results, get rhod and rhop for two lowest bands
    for li, lut in enumerate(luts):

        ## select model based on minimum rmsd between two best fitting bands
        if setu['dsf_model_selection'] == 'min_drmsd':
            if verbosity > 1: print('Computing RMSD for model {}'.format(lut))

            rhop_f = np.zeros((aot_stack[lut]['b1'].shape[0],aot_stack[lut]['b1'].shape[1],2)) + np.nan
            rhod_f = np.zeros((aot_stack[lut]['b1'].shape[0],aot_stack[lut]['b1'].shape[1],2)) + np.nan
            for bi, b in enumerate(aot_bands):
                aot_sub = np.where(aot_stack[lut]['b1']==bi)
                ## get rhod for b1
                if (setu['dsf_path_reflectance'] == 'resolved'):
                    rhod_f[aot_sub[0], aot_sub[1], 0] = gem['data'][gem['bands'][b]['rhot_ds']][aot_sub]
                else:
                    rhod_f[aot_sub[0], aot_sub[1], 0] = dsf_rhod[b][aot_sub]
                ## get rho path for b1
                if len(aot_sub[0]) > 0:
                    if (use_revlut):
                        xi = [gem['data']['pressure'+gk][aot_sub],
                                          gem['data']['raa'+gk][aot_sub],
                                          gem['data']['vza'+gk][aot_sub],
                                          gem['data']['sza'+gk][aot_sub],
                                          gem['data']['wind'+gk][aot_sub]]
                    else:
                        xi = [gem['data']['pressure'+gk],
                                          gem['data']['raa'+gk],
                                          gem['data']['vza'+gk],
                                          gem['data']['sza'+gk],
                                          gem['data']['wind'+gk]]
                    rhop_f[aot_sub[0], aot_sub[1], 0] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd'][par],
                                                                xi[1], xi[2], xi[3], xi[4], aot_stack[lut]['min'][aot_sub]))

                ## get rhod for b1
                aot_sub = np.where(aot_stack[lut]['b2']==bi)
                if (setu['dsf_path_reflectance'] == 'resolved'):
                    rhod_f[aot_sub[0], aot_sub[1], 1] = gem['data'][gem['bands'][b]['rhot_ds']][aot_sub]
                else:
                    rhod_f[aot_sub[0], aot_sub[1], 1] = dsf_rhod[b][aot_sub]
                ## get rhop for b1
                if len(aot_sub[0]) > 0:
                    if (use_revlut):
                        xi = [gem['data']['pressure'+gk][aot_sub],
                                          gem['data']['raa'+gk][aot_sub],
                                          gem['data']['vza'+gk][aot_sub],
                                          gem['data']['sza'+gk][aot_sub],
                                          gem['data']['wind'+gk][aot_sub]]
                    else:
                        xi = [gem['data']['pressure'+gk],
                                          gem['data']['raa'+gk],
                                          gem['data']['vza'+gk],
                                          gem['data']['sza'+gk],
                                          gem['data']['wind'+gk]]
                    rhop_f[aot_sub[0], aot_sub[1], 1] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd'][par],
                                                                xi[1], xi[2], xi[3], xi[4], aot_stack[lut]['min'][aot_sub]))

            ## rmsd for current bands
            cur_sel_par = np.sqrt(np.nanmean(np.square((rhod_f-rhop_f)), axis=2))
        ## end select with min RMSD

        ## select model based on minimum delta tau between two lowest aot bands
        if setu['dsf_model_selection'] == 'min_dtau':
            cur_sel_par = aot_stack[lut]['dtau']
        ## end select with min delta tau

        ## store minimum info
        if li == 0:
            aot_lut = np.zeros(aot_stack[lut]['min'].shape).astype(int)
            aot_lut[aot_stack[lut]['mask']] = -1
            aot_sel = aot_stack[lut]['min'] * 1.0
            aot_sel_par = cur_sel_par * 1.0
        else:
            aot_sub = np.where(cur_sel_par<aot_sel_par)
            if len(aot_sub[0]) == 0: continue
            aot_lut[aot_sub] = li
            aot_sel[aot_sub] = aot_stack[lut]['min'][aot_sub]*1.0
            aot_sel_par[aot_sub] = cur_sel_par[aot_sub] * 1.0

    rhod_f = None
    rhod_p = None

    gem['data']['aot_550'] = aot_sel

    ## set up interpolator for tiled processing
    if setu['dsf_path_reflectance'] == 'tiled':
        xnew = np.linspace(0, tiles[-1][1], gem['gatts']['data_dimensions'][1])
        ynew = np.linspace(0, tiles[-1][0], gem['gatts']['data_dimensions'][0])

    ## compute surface reflectances
    for bi, b in enumerate(gem['bands']):
        if ('rhot_ds' not in gem['bands'][b]) or ('tt_gas' not in gem['bands'][b]): continue
        if gem['bands'][b]['tt_gas'] < setu['min_tgas_rho']: continue
        if gem['bands'][b]['rhot_ds'] not in gem['data']: continue

        dsi = gem['bands'][b]['rhot_ds']
        dso = gem['bands'][b]['rhos_ds']

        t0 = time.time()
        if verbosity > 1: print('Computing surface reflectance', b, gem['bands'][b]['wave_name'], '{:.3f}'.format(gem['bands'][b]['tt_gas']))
        gem['data'][dso] = np.zeros(gem['data'][dsi].shape)+np.nan

        romix = np.zeros(gem['data']['aot_550'].shape)+np.nan
        astot = np.zeros(gem['data']['aot_550'].shape)+np.nan
        dutott = np.zeros(gem['data']['aot_550'].shape)+np.nan

        for li, lut in enumerate(luts):
            ls = np.where(aot_lut == li)
            if len(ls[0]) == 0: continue
            ai = aot_sel[ls]

            ## take all pixels if using fixed processing
            #if aot_lut.shape == (1,1): ls = np.where(gem['data'][dsi])

            if (use_revlut):
                xi = [gem['data']['pressure'+gk][ls],
                      gem['data']['raa'+gk][ls],
                      gem['data']['vza'+gk][ls],
                      gem['data']['sza'+gk][ls],
                      gem['data']['wind'+gk][ls]]
            else:
                xi = [gem['data']['pressure'+gk],
                      gem['data']['raa'+gk],
                      gem['data']['vza'+gk],
                      gem['data']['sza'+gk],
                      gem['data']['wind'+gk]]

            ## path reflectance
            romix[ls] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd'][par], xi[1], xi[2], xi[3], xi[4], ai))
            ## transmittance and spherical albedo
            astot[ls] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd']['astot'], xi[1], xi[2], xi[3], xi[4], ai))
            dutott[ls] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd']['dutott'], xi[1], xi[2], xi[3], xi[4], ai))

        ## interpolate tiled processing to full scene
        if setu['dsf_path_reflectance'] == 'tiled':
            if verbosity > 1: print('Interpolating tiles')
            romix = ac.shared.tiles_interp(romix, xnew, ynew, target_mask=None, smooth=True, kern_size=3, method='linear')
            astot = ac.shared.tiles_interp(astot, xnew, ynew, target_mask=None, smooth=True, kern_size=3, method='linear')
            dutott = ac.shared.tiles_interp(dutott, xnew, ynew, target_mask=None, smooth=True, kern_size=3, method='linear')

        rhot_noatm = (gem['data'][dsi]/ gem['bands'][b]['tt_gas']) - romix
        romix = None
        gem['data'][dso] = (rhot_noatm) / (dutott + astot*rhot_noatm)
        astot=None
        dutott=None
        rhot_noatm = None
        if verbosity > 1: print('{}/B{} took {:.1f}s ({})'.format(gem['gatts']['sensor'], b, time.time()-t0, 'RevLUT' if use_revlut else 'StdLUT'))

    ## compute l8 orange band
    l8_orange_band = True
    if (gem['gatts']['sensor'] == 'L8_OLI') & (l8_orange_band):
        if verbosity > 1: print('Computing orange band')
        ## load orange band configuration
        ob_cfg = ac.shared.import_config(ac.config['data_dir']+'/L8/oli_orange.cfg')

        ## read rsr for wavelength name
        sensor_o = 'L8_OLI_ORANGE'
        rsrd_o = ac.shared.rsr_dict(sensor_o)[sensor_o]
        ob = {k:rsrd_o[k]['O'] for k in ['wave_mu', 'wave_nm', 'wave_name']}
        ob['rhos_ds'] = 'rhos_{}'.format(ob['wave_name'])
        gem['bands']['O'] = ob

        ## compute orange band
        ob_data = gem['data'][gem['bands']['8']['rhos_ds']]*float(ob_cfg['pf'])
        ob_data += gem['data'][gem['bands']['3']['rhos_ds']]*float(ob_cfg['gf'])
        ob_data += gem['data'][gem['bands']['4']['rhos_ds']]*float(ob_cfg['rf'])
        gem['data'][ob['rhos_ds']] = ob_data
        ob_data = None
        ob = None

    ofile = None
    if output_file:
        new_nc = True
        if target_file is None:
            ofile = gemf.replace('_L1R.nc', '_L2R.nc')
            if ('output' in setu) & (output is None): output = setu['output']
            if output is not None: ofile = '{}/{}'.format(output, os.path.basename(ofile))
        else:
            ofile = '{}'.format(target_file)

        ## add settings to gatts
        for k in setu:
            if k in gem['gatts']: continue
            if setu[k] in [True, False]:
                gem['gatts'][k] = str(setu[k])
            else:
                gem['gatts'][k] = setu[k]

        ## write output data
        for bi, b in enumerate(gem['bands']):
            dso = gem['bands'][b]['rhos_ds']
            if dso not in gem['data']: continue
            if verbosity > 1: print('Writing B{} {}'.format(b, dso))
            ds_att = gem['bands'][b]
            ds_att['wavelength']=ds_att['wave_nm']
            ac.output.nc_write(ofile,
                               dso,
                               gem['data'][dso],
                               attributes = gem['gatts'],
                               dataset_attributes = ds_att, new=new_nc)
            new_nc = False

        ## reformat & save aot
        if setu['dsf_path_reflectance'] == 'fixed':
            gem['data']['aot_550'] = np.repeat(gem['data']['aot_550'], gem['gatts']['data_elements']).reshape(gem['gatts']['data_dimensions'])
        if setu['dsf_path_reflectance'] == 'tiled':
            gem['data']['aot_550'] = ac.shared.tiles_interp(gem['data']['aot_550'], xnew, ynew, target_mask=None, smooth=True, kern_size=3, method='linear')
        ## write aot
        ac.output.nc_write(ofile, 'aot_550', gem['data']['aot_550'], attributes = gem['gatts'], new=new_nc)

        ## copy datasets from inputfile
        copy_datasets = setu['copy_datasets']
        if copy_datasets is not None:
            ## copy rhot all from L1R
            if 'rhot_*' in copy_datasets:
                copy_datasets.remove('rhot_*')
                copy_datasets += [ds for ds in gem['datasets'] if ('rhot_' in ds) & (ds not in copy_datasets)]
                print(copy_datasets)
                print()
            ## copy datasets to L2R
            for ds in copy_datasets:
                if (ds not in gem['datasets']):
                    if verbosity > 2: print('{} not found in {}'.format(ds, gemf))
                    continue
                if verbosity > 1: print('Writing {}'.format(ds))
                if ds not in gem['data']:
                    d, da = ac.shared.nc_data(gemf, ds, attributes=True)
                    ac.output.nc_write(ofile, ds, d.data,
                                                dataset_attributes = da,
                                                attributes = gem['gatts'], new=new_nc)
                else:
                    ac.output.nc_write(ofile, ds, gem['data'][ds],
                                                dataset_attributes = gem['atts'][ds],
                                                attributes = gem['gatts'], new=new_nc)

        if verbosity>0: print('Wrote {}'.format(ofile))

    if return_gem:
        return(gem)
    else:
        return(ofile)
