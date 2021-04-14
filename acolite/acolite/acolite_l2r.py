## def acolite_l2r
## runs ACOLITE/DSF on generic extracted miniscene
## written by Quinten Vanhellemont, RBINS
## 2021-03-01
## modifications: 2021-03-11 (QV) forked from acolite_gem

def acolite_l2r(gem,
                output = None,
                sub = None,
                settings = None,

                output_file = True,
                target_file = None,

                return_gem = False,

                verbosity=0):

    import os, time, datetime
    import numpy as np
    import scipy.ndimage
    import acolite as ac

    ## for gaussian smoothing of aot
    from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
    from astropy.convolution import convolve

    time_start = datetime.datetime.now()

    ## read gem file if NetCDF
    if type(gem) is str: gem = ac.gem.gem(gem)
    gemf = gem.file

    ## combine default and user defined settings
    setu = ac.acolite.settings.parse(gem.gatts['sensor'], settings=settings)
    if 'verbosity' in setu: verbosity = setu['verbosity']
    if 'runid' not in setu: setu['runid'] = time_start.strftime('%Y%m%d_%H%M%S')

    ## check blackfill
    if setu['blackfill_skip']:
        rhot_ds = [ds for ds in gem.datasets if 'rhot_' in ds]
        rhot_wv = [int(ds.split('_')[1]) for ds in rhot_ds]
        bi, bw = ac.shared.closest_idx(rhot_wv, setu['blackfill_wave'])
        band_data = 1.0*gem.data(rhot_ds[bi])
        npx = band_data.shape[0] * band_data.shape[1]
        #nbf = npx - len(np.where(np.isfinite(band_data))[0])
        nbf = npx - len(np.where(np.isfinite(band_data)*(band_data>0))[0])
        band_data = None
        if (nbf/npx) >= float(setu['blackfill_max']):
            if verbosity>1: print('Skipping scene as crop is {:.0f}% blackfill'.format(100*nbf/npx))
            return()

    if verbosity > 0: print('Running acolite for {}'.format(gemf))

    output_name = gem.gatts['output_name'] if 'output_name' in gem.gatts else os.path.basename(gemf).replace('.nc', '')

    ## get dimensions and number of elements
    gem.gatts['data_dimensions'] = gem.data('lon').shape
    gem.gatts['data_elements'] = gem.gatts['data_dimensions'][0]*gem.gatts['data_dimensions'][1]

    ## read rsrd and get band wavelengths
    rsrd = ac.shared.rsr_dict(gem.gatts['sensor'])
    if gem.gatts['sensor'] in rsrd:
        rsrd = rsrd[gem.gatts['sensor']]
    else:
        rsrd = None
        stop

    ## set defaults
    gem.gatts['uoz'] = setu['uoz_default']
    gem.gatts['uwv'] = setu['uwv_default']
    gem.gatts['wind'] = setu['wind']
    gem.gatts['pressure'] = setu['pressure']

    ## read ancillary data
    if setu['ancillary_data']:
        clon = np.nanmedian(gem.data('lon'))
        clat = np.nanmedian(gem.data('lat'))
        anc = ac.ac.ancillary.get(gem.gatts['isodate'], clon, clat)

        ## overwrite the defaults
        if ('ozone' in anc): gem.gatts['uoz'] = anc['ozone']['interp']/1000. ## convert from MET data
        if ('p_water' in anc): gem.gatts['uwv'] = anc['p_water']['interp']/10. ## convert from MET data
        if ('z_wind' in anc) & ('m_wind' in anc):
            gem.gatts['wind'] = ((anc['z_wind']['interp']**2) + (anc['m_wind']['interp']**2))**0.5
        if ('press' in anc): gem.gatts['pressure'] = anc['press']['interp']

    ## dem pressure
    if setu['dem_pressure']:
        print('Extracting SRTM DEM data')
        dem = ac.dem.hgt_lonlat(gem.data('lon'), gem.data('lat'))
        dem_pressure = ac.ac.pressure_elevation(dem)
        if setu['dem_pressure_resolved']:
            gem.data_mem['pressure'] = dem_pressure
        else:
            gem.data_mem['pressure'] = np.nanpercentile(dem_pressure, setu['dem_pressure_percentile'])
        gem.datasets.append('pressure')

        if setu['dem_pressure_write']:
            gem.data_mem['dem'] = dem.astype(np.float32)
            gem.data_mem['dem_pressure'] = dem_pressure
        dem = None
        dem_pressure = None

    ## set wind to wind range
    gem.gatts['wind'] = max(0.1, gem.gatts['wind'])
    gem.gatts['wind'] = min(20, gem.gatts['wind'])

    ## get mean average geometry
    geom_ds = ['sza', 'vza', 'raa', 'pressure', 'wind']
    for ds in geom_ds: gem.data(ds, store=True, return_data=False)
    geom_mean = {k: np.nanmean(gem.data(k)) if k in gem.datasets else gem.gatts[k] for k in geom_ds}

    ## get gas transmittance
    tg_dict = ac.ac.gas_transmittance(geom_mean['sza'], geom_mean['vza'],
                                      uoz=gem.gatts['uoz'], uwv=gem.gatts['uwv'],
                                      sensor=gem.gatts['sensor'])

    ## make bands dataset
    gem.bands = {}
    for bi, b in enumerate(rsrd['rsr_bands']):
        if b not in gem.bands:
            gem.bands[b] = {k:rsrd[k][b] for k in ['wave_mu', 'wave_nm', 'wave_name'] if b in rsrd[k]}
            gem.bands[b]['rhot_ds'] = 'rhot_{}'.format(gem.bands[b]['wave_name'])
            gem.bands[b]['rhos_ds'] = 'rhos_{}'.format(gem.bands[b]['wave_name'])
            for k in tg_dict:
                if k not in ['wave']: gem.bands[b][k] = tg_dict[k][b]
            gem.bands[b]['wavelength']=gem.bands[b]['wave_nm']
    ## end bands dataset

    ## atmospheric correction
    if setu['aerosol_correction'] == 'dark_spectrum':
        ac_opt = 'dsf'
    elif  setu['aerosol_correction'] == 'exponential':
        ac_opt = 'exp'
    else:
        print('Option {} {} not configured'.format('aerosol_correction', setu['aerosol_correction']))
        ac_opt = 'dsf'
    print('Using {} atmospheric correction'.format(ac_opt.upper()))

    ## determine use of reverse lut rhot->aot
    use_revlut = False
    ## if path reflectance is tiled or resolved, use reverse lut
    #if setu['dsf_path_reflectance'] != 'fixed': use_revlut = True
    ## no need to use reverse lut if fixed geometry is used
    #if setu['resolved_geometry']:
    ## we want to use the reverse lut to derive aot if the geometry data is resolved
    for ds in geom_ds:
        if ds not in gem.datasets:
            gem.data_mem[ds] = geom_mean[ds]
        else:
            tmp = gem.data(ds, store=True)
        if len(np.atleast_1d(gem.data(ds)))>1:
            use_revlut=True ## if any dataset more than 1 dimension use revlut
        gem.data_mem['{}_mean'.format(ds)] = np.asarray(np.nanmean(gem.data(ds))) ## also store tile mean
        gem.data_mem['{}_mean'.format(ds)].shape+=(1,1) ## make 1,1 dimensions
    if not setu['resolved_geometry']: use_revlut = False

    ## for ease of subsetting later, repeat single element datasets to the tile shape
    if use_revlut:
        for ds in geom_ds:
            if len(np.atleast_1d(gem.data(ds)))!=1: continue
            gem.data_mem[ds] = np.repeat(gem.data_mem[ds], gem.gatts['data_elements']).reshape(gem.gatts['data_dimensions'])
    else:
        ## if use revlut is False at this point, we don't have resolved geometry
        setu['resolved_geometry'] = False
    #print(use_revlut, setu['dsf_path_reflectance'], setu['resolved_geometry'])
    ## end determine revlut

    ## for tiled processing track tile positions and average geometry
    tiles = []
    if 'dsf_tile_dimensions' not in setu: setu['dsf_tile_dimensions'] = None
    if (setu['dsf_path_reflectance'] == 'tiled') & (setu['dsf_tile_dimensions'] is not None):
        ni = np.ceil(gem.gatts['data_dimensions'][0]/setu['dsf_tile_dimensions'][0]).astype(int)
        nj = np.ceil(gem.gatts['data_dimensions'][1]/setu['dsf_tile_dimensions'][1]).astype(int)
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
                    subti[1] = np.min((subti[1], gem.gatts['data_dimensions'][0]))
                    subtj = [setu['dsf_tile_dimensions'][1]*tj, setu['dsf_tile_dimensions'][1]*(tj+1)]
                    subtj[1] = np.min((subtj[1], gem.gatts['data_dimensions'][1]))
                    tiles.append((ti, tj, subti, subtj))

            ## create tile geometry datasets
            for ds in geom_ds:
                if len(np.atleast_1d(gem.data(ds)))>1: ## if not fixed geometry
                    gem.data_mem['{}_tiled'.format(ds)] = np.zeros((ni,nj))+np.nan
                    for t in range(ntiles):
                        ti, tj, subti, subtj = tiles[t]
                        gem.data_mem['{}_tiled'.format(ds)][ti, tj] = \
                            np.nanmean(gem.data(ds)[subti[0]:subti[1],subtj[0]:subtj[1]])
                else: ## if fixed geometry
                    gem.data_mem['{}_tiled'.format(ds)] = gem.data(ds)
    ## end tiling

    ## read LUTs
    if setu['sky_correction']:
        par = 'romix+rsky_t'
    else:
        par = 'romix'

    ## set output directory
    if setu['output'] is not None:
        output_ = setu['output']
    else:
        output_ = os.path.dirname(gemf)

    ## write settings
    settings_file = '{}/acolite_run_{}_l2r_settings.txt'.format(output_,setu['runid'])
    ac.acolite.settings.write(settings_file, setu)
    print(settings_file)

    ## setup output file
    ofile = None
    if output_file:
        new_nc = True
        if target_file is None:
            ofile = gemf.replace('_L1R.nc', '_L2R.nc')
            #if ('output' in setu) & (output is None): output = setu['output']
            if output is None: output = output_
            ofile = '{}/{}'.format(output, os.path.basename(ofile))
        else:
            ofile = '{}'.format(target_file)

        gemo = ac.gem.gem(ofile, new=True)
        gemo.bands = gem.bands
        gemo.verbosity = setu['verbosity']
        gemo.gatts = {k: gem.gatts[k] for k in gem.gatts}
        ## add settings to gatts
        for k in setu:
            if k in gem.gatts: continue
            if setu[k] in [True, False]:
                gemo.gatts[k] = str(setu[k])
            else:
                gemo.gatts[k] = setu[k]

        ## copy datasets from inputfile
        copy_rhot = False
        copy_datasets = setu['copy_datasets']
        if copy_datasets is not None:
            ## copy rhot all from L1R
            if 'rhot_*' in copy_datasets:
                copy_datasets.remove('rhot_*')
                copy_rhot = True
                #copy_datasets += [ds for ds in gem.datasets if ('rhot_' in ds) & (ds not in copy_datasets)]
            ## copy datasets to L2R
            for ds in copy_datasets:
                if (ds not in gem.datasets):
                    if verbosity > 2: print('{} not found in {}'.format(ds, gemf))
                    continue
                if verbosity > 1: print('Writing {}'.format(ds))
                cdata, catts = gem.data(ds, attributes=True)
                gemo.write(ds, cdata, ds_att=catts)

        ## write dem
        if setu['dem_pressure_write']:
            for k in ['dem', 'dem_pressure']:
                if k in gem.data_mem:
                    gemo.write(k, gem.data_mem[k])
                    gem.data_mem[k] = None

    ## load reverse lut romix -> aot
    if use_revlut: revl = ac.aerlut.reverse_lut(gem.gatts['sensor'], par=par)
    ## load aot -> atmospheric parameters lut
    lutdw = ac.aerlut.import_luts(add_rsky=True, sensor=gem.gatts['sensor'])
    luts = list(lutdw.keys())

    ## run through bands to get aot
    aot_dict = {}
    dsf_rhod = {}
    for bi, b in enumerate(gem.bands):
        if (b in setu['dsf_exclude_bands']): continue
        if ('rhot_ds' not in gem.bands[b]) or ('tt_gas' not in gem.bands[b]): continue
        if gem.bands[b]['rhot_ds'] not in gem.datasets: continue

        ## skip band for aot computation
        if gem.bands[b]['tt_gas'] < setu['min_tgas_aot']: continue

        ## skip bands according to configuration
        if (gem.bands[b]['wave_nm'] < setu['dsf_wave_range'][0]): continue
        if (gem.bands[b]['wave_nm'] > setu['dsf_wave_range'][1]): continue
        if (b in setu['dsf_exclude_bands']): continue

        if verbosity > 1: print(b, gem.bands[b]['rhot_ds'])

        band_data = gem.data(gem.bands[b]['rhot_ds'])*1.0
        band_shape = band_data.shape
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
            gk='_mean'
            #if not use_revlut:
            #    gk='_mean'
            #else:
            #    band_data = np.tile(band_data, band_shape)
            if verbosity > 2: print(b, setu['dsf_spectrum_option'], '{:.3f}'.format(band_data[0,0]))

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
        band_data[band_sub] /= gem.bands[b]['tt_gas']

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
                aot_band[lut][band_sub] = revl[lut]['rgi'][b]((gem.data_mem['pressure'+gk][band_sub],
                                                               gem.data_mem['raa'+gk][band_sub],
                                                               gem.data_mem['vza'+gk][band_sub],
                                                               gem.data_mem['sza'+gk][band_sub],
                                                               gem.data_mem['wind'+gk][band_sub],
                                                               band_data[band_sub]))
                # mask out of range aot
                aot_band[lut][aot_band[lut]<=revl[lut]['minaot']]=np.nan
                aot_band[lut][aot_band[lut]>=revl[lut]['maxaot']]=np.nan

            ## standard lut interpolates rhot to results for different aot values
            else:
                ## get rho path for lut steps in aot
                tmp = lutdw[lut]['rgi'][b]((gem.data_mem['pressure'+gk],
                                            lutdw[lut]['ipd'][par],
                                            gem.data_mem['raa'+gk],
                                            gem.data_mem['vza'+gk],
                                            gem.data_mem['sza'+gk],
                                            gem.data_mem['wind'+gk], lutdw[lut]['meta']['tau']))
                tmp = tmp.flatten()

                ## interpolate rho path to observation
                aot_band[lut][band_sub] = np.interp(band_data[band_sub], tmp,
                                                   lutdw[lut]['meta']['tau'],
                                                   left=np.nan, right=np.nan)
                #print(aot_band[lut][band_sub])
            tel = time.time()-t0

            if verbosity > 1: print('{}/B{} {} took {:.3f}s ({})'.format(gem.gatts['sensor'], b, lut, tel, 'RevLUT' if use_revlut else 'StdLUT'))

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

        ## if minimum for fixed retrieval is nan, set it to 0.01
        if setu['dsf_path_reflectance'] == 'fixed':
            if np.isnan(aot_stack[lut]['min']):
                aot_stack[lut]['min'][0][0] = 0.01

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
                ## run through two best fitting bands
                for ai, ab in enumerate(['b1', 'b2']):
                    aot_sub = np.where(aot_stack[lut][ab]==bi)
                    ## get rhod for current band
                    if (setu['dsf_path_reflectance'] == 'resolved'):
                        rhod_f[aot_sub[0], aot_sub[1], ai] = gem.data(gem.bands[b]['rhot_ds'])[aot_sub]
                    else:
                        rhod_f[aot_sub[0], aot_sub[1], ai] = dsf_rhod[b][aot_sub]
                    ## get rho path for current band
                    if len(aot_sub[0]) > 0:
                        if (use_revlut):
                            xi = [gem.data_mem['pressure'+gk][aot_sub],
                                              gem.data_mem['raa'+gk][aot_sub],
                                              gem.data_mem['vza'+gk][aot_sub],
                                              gem.data_mem['sza'+gk][aot_sub],
                                              gem.data_mem['wind'+gk][aot_sub]]
                        else:
                            xi = [gem.data_mem['pressure'+gk],
                                              gem.data_mem['raa'+gk],
                                              gem.data_mem['vza'+gk],
                                              gem.data_mem['sza'+gk],
                                              gem.data_mem['wind'+gk]]
                        rhop_f[aot_sub[0], aot_sub[1], ai] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd'][par],
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

    ## set up interpolator for tiled processing
    if setu['dsf_path_reflectance'] == 'tiled':
        xnew = np.linspace(0, tiles[-1][1], gem.gatts['data_dimensions'][1])
        ynew = np.linspace(0, tiles[-1][0], gem.gatts['data_dimensions'][0])

    ## write aot to outputfile
    if output_file:
        ## reformat & save aot
        if setu['dsf_path_reflectance'] == 'fixed':
            aot_out = np.repeat(aot_sel, gem.gatts['data_elements']).reshape(gem.gatts['data_dimensions'])
        elif setu['dsf_path_reflectance'] == 'tiled':
            aot_out = ac.shared.tiles_interp(aot_sel, xnew, ynew, target_mask=None, smooth=True, kern_size=3, method='linear')
        else:
            aot_out = aot_sel * 1.0
        ## write aot
        gemo.write('aot_550', aot_out)
        aot_out = None

    ## store ttot for glint correction
    ttot_all = {}

    ### store scene mask
    #scene_mask = np.zeros(gemo.gatts['data_dimensions'], dtype=np.uint8)

    print('use_revlut', use_revlut)
    ## compute surface reflectances
    for bi, b in enumerate(gem.bands):
        if ('rhot_ds' not in gem.bands[b]) or ('tt_gas' not in gem.bands[b]): continue
        if gem.bands[b]['rhot_ds'] not in gem.datasets: continue ## skip if we don't have rhot for a band that is in the RSR file

        dsi = gem.bands[b]['rhot_ds']
        dso = gem.bands[b]['rhos_ds']
        cur_data, cur_att = gem.data(dsi, attributes=True)

        ## store rhot in output file
        if copy_rhot:
            gemo.write(dsi, cur_data, ds_att = cur_att)

        if gem.bands[b]['tt_gas'] < setu['min_tgas_rho']: continue
        if gem.bands[b]['rhot_ds'] not in gem.datasets: continue

        t0 = time.time()
        if verbosity > 1: print('Computing surface reflectance', b, gem.bands[b]['wave_name'], '{:.3f}'.format(gem.bands[b]['tt_gas']))

        gem.data_mem[dso] = np.zeros(cur_data.shape)+np.nan
        if setu['slicing']: valid_mask = np.isfinite(cur_data)

        ## shape of atmospheric datasets
        atm_shape = aot_sel.shape
        ## if path reflectance is resolved, but resolved geometry available
        if (use_revlut) & (setu['dsf_path_reflectance'] == 'fixed'):
            atm_shape = cur_data.shape
            gk = ''
        romix = np.zeros(atm_shape)+np.nan
        astot = np.zeros(atm_shape)+np.nan
        dutott = np.zeros(atm_shape)+np.nan
        if setu['glint_correction']:
            ttot_all[b] = np.zeros(atm_shape)+np.nan

        for li, lut in enumerate(luts):
            ls = np.where(aot_lut == li)
            if len(ls[0]) == 0: continue
            ai = aot_sel[ls]

            ## resolved geometry with fixed path reflectance
            if (use_revlut) & (setu['dsf_path_reflectance'] == 'fixed'):
                ls = np.where(cur_data)
            ## take all pixels if using fixed processing
            #if aot_lut.shape == (1,1): ls = np.where(gem['data'][dsi])

            if (use_revlut):
                xi = [gem.data_mem['pressure'+gk][ls],
                      gem.data_mem['raa'+gk][ls],
                      gem.data_mem['vza'+gk][ls],
                      gem.data_mem['sza'+gk][ls],
                      gem.data_mem['wind'+gk][ls]]
            else:
                xi = [gem.data_mem['pressure'+gk],
                      gem.data_mem['raa'+gk],
                      gem.data_mem['vza'+gk],
                      gem.data_mem['sza'+gk],
                      gem.data_mem['wind'+gk]]

            ## path reflectance
            romix[ls] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd'][par], xi[1], xi[2], xi[3], xi[4], ai))
            ## transmittance and spherical albedo
            astot[ls] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd']['astot'], xi[1], xi[2], xi[3], xi[4], ai))
            dutott[ls] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd']['dutott'], xi[1], xi[2], xi[3], xi[4], ai))
            ## total transmittance
            if setu['glint_correction']:
                ttot_all[b][ls] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd']['ttot'], xi[1], xi[2], xi[3], xi[4], ai))

        ## interpolate tiled processing to full scene
        if setu['dsf_path_reflectance'] == 'tiled':
            if verbosity > 1: print('Interpolating tiles')
            romix = ac.shared.tiles_interp(romix, xnew, ynew, target_mask=(valid_mask if setu['slicing'] else None), \
            target_mask_full=True, smooth=True, kern_size=3, method='linear')
            astot = ac.shared.tiles_interp(astot, xnew, ynew, target_mask=(valid_mask if setu['slicing'] else None), \
            target_mask_full=True, smooth=True, kern_size=3, method='linear')
            dutott = ac.shared.tiles_interp(dutott, xnew, ynew, target_mask=(valid_mask if setu['slicing'] else None), \
            target_mask_full=True, smooth=True, kern_size=3, method='linear')
            if setu['glint_correction']:
                ttot_all[b] = ac.shared.tiles_interp(ttot_all[b], xnew, ynew, target_mask=(valid_mask if setu['slicing'] else None), \
                target_mask_full=True, smooth=True, kern_size=3, method='linear')

        ## write ac parameters
        if setu['dsf_write_tiled_parameters']:
            if len(np.atleast_1d(romix)>1):
                if romix.shape == cur_data.shape:
                    gemo.write('romix_{}'.format(gem.bands[b]['wave_name']), romix)
            if len(np.atleast_1d(astot)>1):
                if astot.shape == cur_data.shape:
                    gemo.write('astot_{}'.format(gem.bands[b]['wave_name']), astot)
            if len(np.atleast_1d(dutott)>1):
                if dutott.shape == cur_data.shape:
                    gemo.write('dutott_{}'.format(gem.bands[b]['wave_name']), dutott)
            if setu['glint_correction']:
                if len(np.atleast_1d(ttot_all[b])>1):
                    if ttot_all[b].shape == cur_data.shape:
                        gemo.write('ttot_{}'.format(gem.bands[b]['wave_name']), ttot_all[b])

        ## do atmospheric correction
        rhot_noatm = (cur_data/ gem.bands[b]['tt_gas']) - romix
        romix = None
        cur_data = (rhot_noatm) / (dutott + astot*rhot_noatm)
        astot=None
        dutott=None
        rhot_noatm = None

        ## write rhos
        ds_att = gem.bands[b]
        ds_att['wavelength']=ds_att['wave_nm']
        gemo.write(dso, cur_data, ds_att = ds_att)
        cur_data = None
        if verbosity > 1: print('{}/B{} took {:.1f}s ({})'.format(gem.gatts['sensor'], b, time.time()-t0, 'RevLUT' if use_revlut else 'StdLUT'))
    aot_lut, aot_sel = None, None

    ## update outputfile dataset info
    gemo.datasets_read()

    ## glint correction
    if (setu['aerosol_correction'] == 'dark_spectrum') & setu['glint_correction']:
        ## update output gem datasets to get rhos
        gemo.datasets_read()

        ## find bands for glint correction
        gc_swir1, gc_swir2 = None, None
        gc_swir1_b, gc_swir2_b = None, None
        swir1d, swir2d = 1000, 1000
        gc_user, gc_mask = None, None
        gc_user_b, gc_mask_b = None, None
        userd, maskd = 1000, 1000
        for b in gemo.bands:
            ## swir1
            sd = np.abs(gemo.bands[b]['wave_nm'] - 1600)
            if sd < 100:
                if sd < swir1d:
                    gc_swir1 = gemo.bands[b]['rhos_ds']
                    swir1d = sd
                    gc_swir1_b = b
            ## swir2
            sd = np.abs(gemo.bands[b]['wave_nm'] - 2200)
            if sd < 100:
                if sd < swir2d:
                    gc_swir2 = gemo.bands[b]['rhos_ds']
                    swir2d = sd
                    gc_swir2_b = b
            ## mask band
            sd = np.abs(gemo.bands[b]['wave_nm'] - setu['glint_mask_rhos_wave'])
            if sd < 100:
                if sd < maskd:
                    gc_mask = gemo.bands[b]['rhos_ds']
                    maskd = sd
                    gc_mask_b = b
            ## user band
            if setu['glint_force_band'] is not None:
                sd = np.abs(gemo.bands[b]['wave_nm'] - setu['glint_force_band'])
                if sd < 100:
                    if sd < userd:
                        gc_user = gemo.bands[b]['rhos_ds']
                        userd = sd
                        gc_user_b = b

        ## use user selected  band
        if gc_user is not None:
            gc_swir1, gc_swir1_b = None, None
            gc_swir2, gc_swir2_b = None, None

        ## start glint correction
        if ((gc_swir1 is not None) and (gc_swir2 is not None)) or (gc_user is not None):
            t0 = time.time()
            print('Starting glint correction')

            ## compute scattering angle
            dtor = np.pi / 180.
            sza = gem.data_mem['sza'] * dtor
            vza = gem.data_mem['vza'] * dtor
            raa = gem.data_mem['raa'] * dtor

            muv = np.cos(vza)
            mus = np.cos(sza)
            cos2omega = mus*muv + np.sin(sza)*np.sin(vza)*np.cos(raa)
            omega = np.arccos(np.sqrt(cos2omega))
            omega = np.arccos(cos2omega)/2

            ## read and resample refractive index
            refri = ac.ac.refri()
            refri_sen = ac.shared.rsr_convolute_dict(refri['wave']/1000, refri['n'], rsrd['rsr'])

            ## compute fresnel reflectance for each n
            Rf_sen = {b: ac.ac.sky_refl(omega, n_w=refri_sen[b]) for b in refri_sen}

            ## compute where to apply the glint correction
            ## sub_gc has the idx for non masked data with rhos_ref below the masking threshold
            gc_mask_data = gemo.data(gc_mask)
            sub_gc = np.where(np.isfinite(gc_mask_data) & \
                              (gc_mask_data<=setu['glint_mask_rhos_threshold']))
            gc_mask_data = None

            ## get reference bands transmittance
            for ib, b in enumerate(gemo.bands):
                rhos_ds = gemo.bands[b]['rhos_ds']
                if rhos_ds not in [gc_swir1, gc_swir2, gc_user]: continue
                if rhos_ds not in gemo.datasets: continue

                ## two way direct transmittance
                T_cur  = np.exp(-1.*(ttot_all[b]/muv)) * np.exp(-1.*(ttot_all[b]/mus))
                if rhos_ds == gc_user:
                    T_USER = T_cur[sub_gc]
                else:
                    if rhos_ds == gc_swir1:
                        T_SWIR1 = T_cur[sub_gc]
                    if rhos_ds == gc_swir2:
                        T_SWIR2 = T_cur[sub_gc]
                T_cur = None

            ## swir band choice is made for first band
            gc_choice = False
            ## glint correction per band
            for ib, b in enumerate(gemo.bands):
                rhos_ds = gemo.bands[b]['rhos_ds']
                if rhos_ds not in gemo.datasets: continue
                if b not in ttot_all: continue
                print('Performing glint correction for band {} ({} nm)'.format(b, gemo.bands[b]['wave_name']))

                ## two way direct transmittance
                T_cur  = np.exp(-1.*(ttot_all[b]/muv)) * np.exp(-1.*(ttot_all[b]/mus))

                ## get gc factors for this band
                if gc_user is None:
                    gc_SWIR1 = (T_cur[sub_gc]/T_SWIR1) * (Rf_sen[b][sub_gc]/Rf_sen[gc_swir1_b][sub_gc])
                    gc_SWIR2 = (T_cur[sub_gc]/T_SWIR2) * (Rf_sen[b][sub_gc]/Rf_sen[gc_swir2_b][sub_gc])
                else:
                    gc_USER = (T_cur[sub_gc]/T_USER) * (Rf_sen[b][sub_gc]/Rf_sen[gc_user_b][sub_gc])

                ## choose glint correction band (based on first band results)
                if gc_choice is False:
                    gc_choice = True
                    if gc_user is None:
                        swir1_rhos = gemo.data(gc_swir1)[sub_gc]
                        swir2_rhos = gemo.data(gc_swir2)[sub_gc]
                        ## set negatives to 0
                        swir1_rhos[swir1_rhos<0] = 0
                        swir2_rhos[swir2_rhos<0] = 0
                        ## estimate glint correction in the blue band
                        g1_blue = gc_SWIR1 * swir1_rhos
                        g2_blue = gc_SWIR2 * swir2_rhos
                        ## use SWIR1 or SWIR2 based glint correction
                        use_swir1 = np.where(g1_blue<g2_blue)
                        g1_blue, g2_blue = None, None
                        rhog_ref = swir2_rhos
                        rhog_ref[use_swir1] = swir1_rhos[use_swir1]
                        swir1_rhos, swir2_rhos = None, None
                        use_swir1 = None
                    else:
                        rhog_ref = gemo.data(gc_user)[sub_gc]
                        ## set negatives to 0
                        rhog_ref[rhog_ref<0] = 0
                    ## write reference glint
                    if setu['glint_write_rhog_ref']:
                        tmp = np.zeros(gemo.gatts['data_dimensions'], dtype=np.float32) + np.nan
                        tmp[sub_gc] = rhog_ref
                        gemo.write('rhog_ref', tmp)
                        tmp = None
                ## end select glint correction band

                ## calculate glint in this band
                if gc_user is None:
                    cur_rhog = gc_SWIR2 * rhog_ref
                    try:
                        cur_rhog[use_swir1] = gc_SWIR1[use_swir1] * rhog_ref[use_swir1]
                    except:
                        cur_rhog[use_swir1] = gc_SWIR1 * rhog_ref[use_swir1]
                else:
                    cur_rhog = gc_USER * rhog_ref

                ## remove glint from rhos
                cur_data = gemo.data(rhos_ds)
                cur_data[sub_gc]-=cur_rhog
                gemo.write(rhos_ds, cur_data, ds_att = gem.bands[b])

                ## write band glint
                if setu['glint_write_rhog_all']:
                    tmp = np.zeros(gemo.gatts['data_dimensions'], dtype=np.float32) + np.nan
                    tmp[sub_gc] = cur_rhog
                    gemo.write('rhog_{}'.format(gemo.bands[b]['wave_name']), tmp, ds_att={'wavelength':gemo.bands[b]['wavelength']})
                    tmp = None
                cur_rhog = None
            Rf_sen = None
            rhog_ref = None
    ## end glint correction

    ## compute l8 orange band
    l8_orange_band = True
    if (gemo.gatts['sensor'] == 'L8_OLI') & (l8_orange_band):
        if verbosity > 1: print('Computing orange band')
        ## load orange band configuration
        ob_cfg = ac.shared.import_config(ac.config['data_dir']+'/L8/oli_orange.cfg')

        ## read rsr for wavelength name
        sensor_o = 'L8_OLI_ORANGE'
        rsrd_o = ac.shared.rsr_dict(sensor_o)[sensor_o]
        ob = {k:rsrd_o[k]['O'] for k in ['wave_mu', 'wave_nm', 'wave_name']}
        ob['rhos_ds'] = 'rhos_{}'.format(ob['wave_name'])
        ob['wavelength']=ds_att['wave_nm']
        gemo.bands['O'] = ob

        ## compute orange band
        ob_data = gemo.data(gemo.bands['8']['rhos_ds'])*float(ob_cfg['pf'])
        ob_data += gemo.data(gemo.bands['3']['rhos_ds'])*float(ob_cfg['gf'])
        ob_data += gemo.data(gemo.bands['4']['rhos_ds'])*float(ob_cfg['rf'])
        gemo.write(ob['rhos_ds'], ob_data, ds_att = ob)
        ob_data = None
        ob = None
    ## end orange band

    if verbosity>0: print('Wrote {}'.format(ofile))

    if return_gem:
        return(gem, setu)
    else:
        return(ofile, setu)
