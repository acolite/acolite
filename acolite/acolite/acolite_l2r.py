## def acolite_l2r
## runs ACOLITE/DSF on generic extracted miniscene
## written by Quinten Vanhellemont, RBINS
## 2021-03-01
## modifications: 2021-03-11 (QV) forked from acolite_gem
##                2021-12-08 (QV) added nc_projection
##                2022-01-01 (QV) added segmented dsf option
##                2022-06-21 (QV) moved orange band to separate function
##                2022-07-15 (QV) added option to select most common model for non-fixed DSF
##                2022-09-24 (QV) removed special case for DESIS
##                2023-03-30 (QV) added option to include band name in rhot/rhos datasets
##                2023-07-12 (QV) removed netcdf_compression settings from gem call
##                2023-10-31 (QV) improved memory management, moved ttot_all interpolation to glint correction
##                2023-12-07 (QV) option to use S2 AUX
##                2024-03-14 (QV) update settings handling
##                2024-04-23 (MB) use ancillary data included in resampled MSI
##                2024-07-22 (QV) include PACE_OCI SWIR RSR
##                2024-10-16 (QV) added RSR versioning support
##                2024-12-17 (QV) changed aerosol_correction to atmospheric_correction_method
##                2025-01-31 (QV) check if lat/lon are present for dem
##                2025-02-03 (QV) use downward gas transmittance for output_ed
##                2025-02-04 (QV) updated settings parsing
##                2025-02-10 (QV) added optimisation option
##                2025-02-11 (QV) switch to settings.merge
##                2025-03-04 (QV) update to hyperspectral RSR
##                2025-03-05 (QV) added optional printouts for bands skipped in DSF
##                2025-03-10 (QV) fix hyperspectral model selection, use setu['verbosity']
##                2025-04-30 (QV) added sensor noise bias correction
##                2025-05-16 (QV) added filtering for sensor noise bias correction
##                2025-07-07 (QV) don't do model selection if only one model is used
##                2025-10-06 (QV) added dsf_aot_option=ancillary and dsf_aot_option=ancillary_fixed

def acolite_l2r(gem,
                output = None,
                sub = None,
                settings = None,

                output_file = True,
                target_file = None,

                return_gem = False):

    import os, time, datetime
    import numpy as np
    import scipy.ndimage, scipy.interpolate, scipy.stats
    import acolite as ac
    import skimage.measure

    time_start = datetime.datetime.now()

    ## read gem file if NetCDF
    if type(gem) is str:
        gem = ac.gem.gem(gem)
    gemf = gem.file

    ## get run/user/sensor settings
    setu = ac.acolite.settings.merge(sensor = gem.gatts['sensor'], settings = settings)

    if 'runid' not in setu: setu['runid'] = time_start.strftime('%Y%m%d_%H%M%S')

    ## convert exclude bands to list
    if setu['dsf_exclude_bands'] != None:
        if type(setu['dsf_exclude_bands']) != list:
            setu['dsf_exclude_bands'] = [setu['dsf_exclude_bands']]
    else:
        setu['dsf_exclude_bands'] = []

    rhot_ds = [ds for ds in gem.datasets if 'rhot_' in ds]
    if len(rhot_ds) == 0:
        print('No TOA reflectance data in {}'.format(gemf))
        return()

    ## check blackfill
    if setu['blackfill_skip']:
        rhot_wv = [int(ds.split('_')[-1]) for ds in rhot_ds] ## use last element of rhot name as wavelength
        bi, bw = ac.shared.closest_idx(rhot_wv, setu['blackfill_wave'])
        band_data = 1.0*gem.data(rhot_ds[bi])
        npx = band_data.shape[0] * band_data.shape[1]
        #nbf = npx - len(np.where(np.isfinite(band_data))[0])
        nbf = npx - len(np.where(np.isfinite(band_data)*(band_data>0))[0])
        del band_data
        if (nbf/npx) >= float(setu['blackfill_max']):
            if setu['verbosity']>0: print('Skipping scene as crop is {:.0f}% blackfill'.format(100*nbf/npx))
            return()

    if setu['verbosity'] > 0: print('Running acolite for {}'.format(gemf))

    ## optimised aot
    if setu['dsf_aot_estimate'].startswith('optim'):
        print('Running AOT optimisation to target spectrum for L2R processing')
        ret = ac.ac.optimise_aot_homogeneous(gem, settings = setu)
        if ret is None:
            print('Error in aot optimisation.')
            return
        else:
            opt_lut, opt_aot = ret
        print('Setting dsf_fixed_aot={:.5f} and dsf_fixed_lut={} to optimisation results'.format(opt_aot, opt_lut))
        setu['dsf_fixed_aot'] = opt_aot
        setu['dsf_fixed_lut'] = opt_lut
    ## end optimised aot

    ## ancillary aot
    if setu['dsf_aot_estimate'].startswith('ancillary'):
        print('Retrieving AOT from ancillary data for L2R processing')
        if setu['dsf_aot_estimate'] == 'ancillary':
            aer_anc = ac.ac.ancillary.get_aer(gem.gatts['isodate'], gem.data('lon'), gem.data('lat'))
        elif setu['dsf_aot_estimate'] == 'ancillary_fixed':
            aer_anc = ac.ac.ancillary.get_aer(gem.gatts['isodate'], np.nanmean(gem.data('lon')), np.nanmean(gem.data('lat')))
        else:
            print('dsf_aot_estimate={} not configured.'.format(setu['dsf_aot_estimate']))
            return

        if 'data' not in aer_anc:
            print('Error in ancillary AOT retrieval.')
            return
        elif aer_anc['data'] == {}:
            print('Error in ancillary AOT retrieval.')
            return
        else:
            aer_ang = aer_anc['data']['TOTANGSTR']['interp']
            aer_aot = aer_anc['data']['TOTEXTTAU']['interp']
            aer_ang_mean = np.nanmean(aer_ang)
            ## 6SV models C(ontinental) 1.08, M(aritime) 0.28,
            anc_ang_threshold = 0.68 ## midway between C and M
            anc_lut = setu['luts'][0]
            if len(setu['luts']) > 1:
                print('Selecting LUT based on ancillary mean angstrom {:.2f} from LUTs {}'.format(aer_ang_mean, setu['luts']))
                print('Angstrom threshold M < {} < C'.format(anc_ang_threshold))
                for lut in setu['luts']:
                    if (aer_ang_mean >= anc_ang_threshold) & (lut[-1] == '1'): anc_lut = '{}'.format(lut)
                    elif (aer_ang_mean < anc_ang_threshold) & (lut[-1] == '2'): anc_lut = '{}'.format(lut)
            anc_aot = np.nanmean(aer_aot) * 1
        print('Setting dsf_fixed_aot={:.3f} (mean) and dsf_fixed_lut={} (mean angstrom={:.2f}) based on ancillary data'.format(anc_aot, anc_lut, aer_ang_mean))
        setu['dsf_fixed_aot'] = anc_aot
        setu['dsf_fixed_lut'] = anc_lut
    ## end ancillary aot

    output_name = gem.gatts['output_name'] if 'output_name' in gem.gatts else os.path.basename(gemf).replace('.nc', '')

    ## get dimensions and number of elements
    gem.gatts['data_dimensions'] = gem.data(gem.datasets[-1]).shape
    gem.gatts['data_elements'] = gem.gatts['data_dimensions'][0]*gem.gatts['data_dimensions'][1]

    ## read rsrd and get band wavelengths
    hyper = False
    sensor_version = None
    sensor_lut = None

    ## determine if sensor in hyperspectral sensors
    ## this determines whether the generic LUT is loaded or a sensor specific one
    if gem.gatts['sensor'] in ac.config['hyper_sensors']: hyper = True

    ## create or load RSR
    if (hyper) & ('band_waves' in gem.gatts) & ('band_widths' in gem.gatts):
        ## make hyperspectral RSR
        rsr = ac.shared.rsr_hyper(gem.gatts['band_waves'],
                                  gem.gatts['band_widths'], step=0.1)
        ## update PACE OCI response with known RSR
        if gem.gatts['sensor'] == 'PACE_OCI':
            ## SWIR
            swir_bands = [282, 283, 284, 285, 286, 287, 288, 289, 290]
            rsrd_swir = ac.shared.rsr_dict('PACE_OCI_SWIR')
            for bi, b in enumerate(swir_bands):
                swir_b = rsrd_swir['PACE_OCI_SWIR']['rsr_bands'][bi]
                #print(bi, b, swir_b, rsrd_swir['PACE_OCI_SWIR']['wave_nm'][swir_b])
                rsr[b]['wave'] = np.asarray(rsrd_swir['PACE_OCI_SWIR']['rsr'][swir_b]['wave'])
                rsr[b]['response'] = np.asarray(rsrd_swir['PACE_OCI_SWIR']['rsr'][swir_b]['response'])
            del rsrd_swir
        ## end update PACE OCI response
        rsrd = ac.shared.rsr_dict(rsrd={gem.gatts['sensor']:{'rsr':rsr}})
        del rsr
    else:
        if setu['rsr_version'] is None:
            rsrd = ac.shared.rsr_dict(gem.gatts['sensor'])
        else: ## use versioned RSR
            sensor_version = '{}_{}'.format(gem.gatts['sensor'], setu['rsr_version'])
            rsrd = ac.shared.rsr_dict(sensor_version)

    if gem.gatts['sensor'] in rsrd:
        sensor_lut = '{}'.format(gem.gatts['sensor'])
        rsrd = rsrd[sensor_lut]
    elif (sensor_version is not None) & (sensor_version in rsrd):
        print('Using RSR {} for sensor {}'.format(sensor_version, gem.gatts['sensor']))
        sensor_lut = '{}'.format(sensor_version)
        rsrd = rsrd[sensor_lut]
    else:
        print('Could not find {} RSR'.format(gem.gatts['sensor']))
        return()

    ## set defaults
    gem.gatts['uoz'] = setu['uoz_default']
    gem.gatts['uwv'] = setu['uwv_default']
    gem.gatts['wind'] = setu['wind']
    gem.gatts['pressure'] = setu['pressure']

    ## read ancillary data
    if (setu['ancillary_data']) & ((('lat' in gem.datasets) & ('lon' in gem.datasets)) | (('lat' in gem.gatts) & ('lon' in gem.gatts))):
        if ('lat' in gem.datasets) & ('lon' in gem.datasets):
            clon = np.nanmedian(gem.data('lon'))
            clat = np.nanmedian(gem.data('lat'))
        else:
            clon = gem.gatts['lon']
            clat = gem.gatts['lat']
        print('Getting ancillary data for {} {:.3f}E {:.3f}N'.format(gem.gatts['isodate'], clon, clat))
        anc = ac.ac.ancillary.get(gem.gatts['isodate'], clon, clat, verbosity=setu['verbosity'])
        for k in ['uoz', 'uwv', 'wind', 'pressure']:
            if (k == 'pressure'):
                if (setu['pressure'] != setu['pressure_default']): continue
            if (k == 'wind') & (setu['wind'] is not None): continue
            if k in anc: gem.gatts[k] = 1.0 * anc[k]
        del clon, clat, anc
    else:
        if (setu['s2_auxiliary_default']) & (gem.gatts['sensor'][0:2] == 'S2') & (gem.gatts['sensor'][4:] == 'MSI'):
            ## get mid point values from AUX ECMWFT
            if 'AUX_ECMWFT_msl_values' in gem.gatts:
                gem.gatts['pressure'] = gem.gatts['AUX_ECMWFT_msl_values'][int(len(gem.gatts['AUX_ECMWFT_msl_values'])/2)]/100 # convert from Pa to hPa
            if 'AUX_ECMWFT_tco3_values' in gem.gatts:
                gem.gatts['uoz'] = gem.gatts['AUX_ECMWFT_tco3_values'][int(len(gem.gatts['AUX_ECMWFT_tco3_values'])/2)] ** -1 * 0.0021415 # convert from kg/m2 to cm-atm
            if 'AUX_ECMWFT_tcwv_values' in gem.gatts:
                gem.gatts['uwv'] = gem.gatts['AUX_ECMWFT_tcwv_values'][int(len(gem.gatts['AUX_ECMWFT_tcwv_values'])/2)]/10 # convert from kg/m2 to g/cm2
            if 'AUX_ECMWFT__0u_values' in gem.gatts and 'AUX_ECMWFT__0v_values' in gem.gatts:  # S2Resampling, msiresampling
                u_wind = gem.gatts['AUX_ECMWFT__0u_values'][int(len(gem.gatts['AUX_ECMWFT__0u_values'])/2)]
                v_wind = gem.gatts['AUX_ECMWFT__0v_values'][int(len(gem.gatts['AUX_ECMWFT__0v_values'])/2)]
                gem.gatts['wind'] = np.sqrt(u_wind * u_wind + v_wind * v_wind)

            ## interpolate to scene centre
            if setu['s2_auxiliary_interpolate']:
                if ('lat' in gem.datasets) & ('lon' in gem.datasets):
                    clon = np.nanmedian(gem.data('lon'))
                    clat = np.nanmedian(gem.data('lat'))
                else:
                    clon = gem.gatts['lon']
                    clat = gem.gatts['lat']
                ## pressure
                aux_par = 'ECMWFT_msl'
                if 'AUX_{}_values'.format(aux_par) in gem.gatts:
                    lndi = scipy.interpolate.LinearNDInterpolator((gem.gatts['AUX_{}_longitudes'.format(aux_par)],
                                                                   gem.gatts['AUX_{}_latitudes'.format(aux_par)]),
                                                                   gem.gatts['AUX_{}_values'.format(aux_par)])
                    gem.gatts['pressure'] = lndi(clon, clat)/100 # convert from Pa to hPa
                    del lndi
                elif 'msl' in gem.datasets:  # S2Resampling, msiresampling
                    lndi = gem.data('msl')
                    gem.gatts['pressure'] = lndi(lndi.len/2, lndi.len/2) / 100  # convert from Pa to hPa
                    del lndi
                ## ozone
                aux_par = 'ECMWFT_tco3'
                if 'AUX_{}_values'.format(aux_par) in gem.gatts:
                    lndi = scipy.interpolate.LinearNDInterpolator((gem.gatts['AUX_{}_longitudes'.format(aux_par)],
                                                                   gem.gatts['AUX_{}_latitudes'.format(aux_par)]),
                                                                   gem.gatts['AUX_{}_values'.format(aux_par)])
                    gem.gatts['uoz'] = lndi(clon, clat)**-1 * 0.0021415 # convert from kg/m2 to cm-atm
                    del lndi
                ## water vapour
                aux_par = 'ECMWFT_tcwv'
                if 'AUX_{}_values'.format(aux_par) in gem.gatts:
                    lndi = scipy.interpolate.LinearNDInterpolator((gem.gatts['AUX_{}_longitudes'.format(aux_par)],
                                                                   gem.gatts['AUX_{}_latitudes'.format(aux_par)]),
                                                                   gem.gatts['AUX_{}_values'.format(aux_par)])
                    gem.gatts['uwv'] = lndi(clon, clat)/10 # convert from kg/m2 to g/cm2
                    del lndi
                elif 'tcwv' in gem.datasets:
                    lndi = gem.data('tcwv')
                    gem.gatts['uwv'] = lndi(lndi.len/2, lndi.len/2) / 10  # convert from kg/m2 to g/cm2
                    del lndi
                del clon, clat

    ## elevation provided
    if setu['elevation'] is not None:
        setu['pressure'] = ac.ac.pressure_elevation(setu['elevation'])
        gem.gatts['pressure'] = setu['pressure']
        print(gem.gatts['pressure'])

    ## dem pressure
    if (setu['dem_pressure']) & ((('lat' in gem.datasets) & ('lon' in gem.datasets)) | (('lat' in gem.gatts) & ('lon' in gem.gatts))):
        if setu['verbosity'] > 1: print('Extracting {} DEM data'.format(setu['dem_source']))
        if ('lat' in gem.datasets) & ('lon' in gem.datasets):
            dem = ac.dem.dem_lonlat(gem.data('lon'), gem.data('lat'), source = setu['dem_source'])
        else:
            dem = ac.dem.dem_lonlat(gem.gatts['lon'], gem.gatts['lon'], source = setu['dem_source'])

        if dem is not None:
            dem_pressure = ac.ac.pressure_elevation(dem)

            if setu['dem_pressure_resolved']:
                gem.data_mem['pressure'] = dem_pressure
                gem.gatts['pressure'] = np.nanmean(gem.data_mem['pressure'])
            else:
                gem.data_mem['pressure'] = np.nanpercentile(dem_pressure, setu['dem_pressure_percentile'])
                gem.gatts['pressure'] = gem.data_mem['pressure']
            gem.datasets.append('pressure')

            if setu['dem_pressure_write']:
                gem.data_mem['dem'] = dem.astype(np.float32)
                gem.data_mem['dem_pressure'] = dem_pressure
        else:
            if setu['verbosity'] > 1: print('Could not determine elevation from {} DEM data'.format(setu['dem_source']))

        dem = None
        dem_pressure = None

    print('default uoz: {:.2f} uwv: {:.2f} pressure: {:.2f}'.format(setu['uoz_default'], setu['uwv_default'], setu['pressure_default']))
    print('current uoz: {:.2f} uwv: {:.2f} pressure: {:.2f}'.format(gem.gatts['uoz'], gem.gatts['uwv'], gem.gatts['pressure']))

    ## which LUT data to read
    par = 'romix'
    ipar = 'romix+rsky_t'
    if (setu['dsf_interface_reflectance']):
        if (setu['dsf_interface_option'] == 'default'):
            par = 'romix+rsky_t'
        elif (setu['dsf_interface_option']  == '6sv'):
            par = 'romix+rsurf'
            ipar = 'romix+rsurf'
        elif (setu['dsf_interface_option']  == 'ffss_boa'):
            ## keep default par/ipar
            ## import skydome
            meta_rsky, lut_rsky, rgi_rsky = ac.ac.skydome.import_skydome_lut(sensor = sensor_lut if not hyper else None, lut_base_skydome = 'ACOLITE-FFSS-202511-82W', par = 'rsky')
        elif (setu['dsf_interface_option']  == 'ffss_toa'):
            par = 'romix+ffss_toa'
            ipar = 'romix+ffss_toa'
        else:
            print('dsf_interface_option={} not configured'.format(setu['dsf_interface_option']))
            return

    ## set wind to wind range
    if gem.gatts['wind'] is None: gem.gatts['wind'] = setu['wind_default']
    if (par == 'romix+rsurf'):
        gem.gatts['wind'] = max(2, gem.gatts['wind'])
        gem.gatts['wind'] = min(20, gem.gatts['wind'])
    elif (par == 'romix+rsky_t') | (par == 'romix'):
        gem.gatts['wind'] = max(0.1, gem.gatts['wind'])
        gem.gatts['wind'] = min(20, gem.gatts['wind'])
    else:
        ## placeholder for ffss_toa
        gem.gatts['wind'] = 0

    ## get mean average geometry
    geom_ds = ['sza', 'vza', 'raa', 'pressure', 'wind']
    for ds in gem.datasets:
        if ('raa_' in ds) or ('vza_' in ds):
            gem.data(ds, store=True, return_data=False)
            geom_ds.append(ds)
    for ds in geom_ds:
        gem.data(ds, store=True, return_data=False)
        if (ds == 'sza'):
            sza = gem.data(ds, store=True, return_data=True)
            if sza is not None:
                high_sza = np.where(sza>setu['sza_limit'])
                if len(high_sza[0]) > 0:
                    print('Warning: SZA out of LUT range')
                    print('Mean SZA: {:.3f}'.format(np.nanmean(sza)))
                    if (setu['sza_limit_replace']):
                        sza[high_sza] = setu['sza_limit']
                        print('Mean SZA after replacing SZA > {}: {:.3f}'.format(setu['sza_limit'],np.nanmean(sza)))
                del sza, high_sza
    geom_mean = {k: np.nanmean(gem.data(k)) if k in gem.datasets else gem.gatts[k] for k in geom_ds}
    if (geom_mean['sza'] > setu['sza_limit']):
        print('Warning: SZA out of LUT range')
        print('Mean SZA: {:.3f}'.format(geom_mean['sza']))
        if  (setu['sza_limit_replace']):
            geom_mean['sza'] = setu['sza_limit']
            print('Mean SZA after replacing SZA > {}: {:.3f}'.format(setu['sza_limit'],geom_mean['sza']))

    if (geom_mean['vza'] > setu['vza_limit']):
        print('Warning: VZA out of LUT range')
        print('Mean VZA: {:.3f}'.format(geom_mean['vza']))
        if  (setu['vza_limit_replace']):
            geom_mean['vza'] = setu['vza_limit']
            print('Mean VZA after replacing VZA > {}: {:.3f}'.format(setu['vza_limit'],geom_mean['vza']))

    ## get gas transmittance
    tg_dict = ac.ac.gas_transmittance(geom_mean['sza'], geom_mean['vza'],
                                      uoz=gem.gatts['uoz'], uwv=gem.gatts['uwv'],
                                      rsr=rsrd['rsr'])

    ## for Ed
    if setu['output_ed']:
        ## solar irradiance
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data']), rsrd['rsr'])
        ## downward gas transmittance
        tdgas = ac.ac.gaslut_interp(geom_mean['sza'], geom_mean['vza'], pressure=geom_mean['pressure'],
                                    sensor=None, lutconfig='202402F', pars = ['dtdica','dtoxyg','dtniox','dtmeth'])
        tdg = np.prod([tdgas[k] for k in tdgas if k != 'wave'], axis=0)
        tdg_b = ac.shared.rsr_convolute_dict(np.asarray(tdgas['wave']), tdg, rsrd['rsr'])

    ## make bands dataset
    gem.bands = {}
    for bi, b in enumerate(rsrd['rsr_bands']):
        if b not in gem.bands:
            gem.bands[b] = {k:rsrd[k][b] for k in ['wave_mu', 'wave_nm', 'wave_name'] if b in rsrd[k]}
            gem.bands[b]['rhot_ds'] = 'rhot_{}'.format(gem.bands[b]['wave_name'])
            gem.bands[b]['rhos_ds'] = 'rhos_{}'.format(gem.bands[b]['wave_name'])
            if setu['add_band_name']:
                gem.bands[b]['rhot_ds'] = 'rhot_{}_{}'.format(b, gem.bands[b]['wave_name'])
                gem.bands[b]['rhos_ds'] = 'rhos_{}_{}'.format(b, gem.bands[b]['wave_name'])
            if setu['add_detector_name']:
                dsname = rhot_ds[bi][5:]
                gem.bands[b]['rhot_ds'] = 'rhot_{}'.format(dsname)
                gem.bands[b]['rhos_ds'] = 'rhos_{}'.format(dsname)
            if setu['output_ed']:
                gem.bands[b]['F0'] = f0_b[b]
                gem.bands[b]['td_gas'] = tdg_b[b]
            for k in tg_dict:
                if k not in ['wave']:
                    gem.bands[b][k] = tg_dict[k][b]
                    if setu['gas_transmittance'] is False: gem.bands[b][k] = 1.0
            gem.bands[b]['wavelength']=gem.bands[b]['wave_nm']

            ## add sensor noise to band attributes
            if (setu['sensor_noise_bias_correction']) & (setu['sensor_noise'] is not None):
                if len(setu['sensor_noise']) == len(rsrd['rsr_bands']):
                    gem.bands[b]['sensor_noise'] = setu['sensor_noise'][bi]
    ## end bands dataset
    del rhot_ds, tg_dict

    ## atmospheric correction
    if setu['atmospheric_correction_method'] == 'dark_spectrum':
        ac_opt = 'dsf'
    elif  setu['atmospheric_correction_method'] == 'exponential':
        ac_opt = 'exp'
    else:
        print('Option {} {} not configured'.format('atmospheric_correction_method', setu['atmospheric_correction_method']))
        ac_opt = 'dsf'
    print('Using {} atmospheric correction'.format(ac_opt.upper()))

    ## determine use of reverse lut rhot->aot
    if setu['resolved_geometry'] & hyper:
        print('Resolved geometry for hyperspectral sensors currently not supported')
        setu['resolved_geometry'] = False

    use_revlut = False
    per_pixel_geometry = False
    ## if path reflectance is tiled or resolved, use reverse lut
    ## no need to use reverse lut if fixed geometry is used
    ## we want to use the reverse lut to derive aot if the geometry data is resolved
    for ds in geom_ds:
        if ds not in gem.datasets:
            gem.data_mem[ds] = geom_mean[ds]
        else:
            tmp = gem.data(ds, store=True)
            del tmp
        if len(np.atleast_1d(gem.data(ds)))>1:
            use_revlut=True ## if any dataset more than 1 dimension use revlut
            per_pixel_geometry = True
        else: ## convert floats into arrays
            gem.data_mem[ds] = np.asarray(gem.data_mem[ds])
            gem.data_mem[ds].shape+=(1,1)
        gem.data_mem['{}_mean'.format(ds)] = np.asarray(np.nanmean(gem.data(ds))) ## also store tile mean
        gem.data_mem['{}_mean'.format(ds)].shape+=(1,1) ## make 1,1 dimensions

    del geom_mean

    ## for tiled processing track tile positions and average geometry
    tiles = []
    if 'dsf_tile_dimensions' not in setu: setu['dsf_tile_dimensions'] = None
    if (setu['dsf_aot_estimate'] == 'tiled') & (setu['dsf_tile_dimensions'] is not None):
        ni = np.ceil(gem.gatts['data_dimensions'][0]/setu['dsf_tile_dimensions'][0]).astype(int)
        nj = np.ceil(gem.gatts['data_dimensions'][1]/setu['dsf_tile_dimensions'][1]).astype(int)
        if (ni <= 1) | (nj <= 1):
            if setu['verbosity'] > 1: print('Scene too small for tiling ({}x{} tiles of {}x{} pixels), using fixed processing'.format(ni,nj,setu['dsf_tile_dimensions'][0], setu['dsf_tile_dimensions'][1]))
            setu['dsf_aot_estimate'] = 'fixed'
        else:
            ntiles = ni*nj
            if setu['verbosity'] > 1: print('Processing with {} tiles ({}x{} tiles of {}x{} pixels)'.format(ntiles, ni, nj, setu['dsf_tile_dimensions'][0], setu['dsf_tile_dimensions'][1]))

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
                    gem.data_mem['{}_tiled'.format(ds)] = np.zeros((ni,nj), dtype=np.float32)+np.nan
                    for t in range(ntiles):
                        ti, tj, subti, subtj = tiles[t]
                        gem.data_mem['{}_tiled'.format(ds)][ti, tj] = \
                            np.nanmean(gem.data(ds)[subti[0]:subti[1],subtj[0]:subtj[1]])
                else: ## if fixed geometry
                    if per_pixel_geometry:
                        gem.data_mem['{}_tiled'.format(ds)] = np.zeros((ni,nj), dtype=np.float32)+gem.data(ds)
                    else:
                        gem.data_mem['{}_tiled'.format(ds)] = 1.0 * gem.data(ds)

            ## remove full geom datasets as tiled will be used
            for ds in geom_ds:
                if '{}_tiled'.format(ds) in gem.data_mem:
                    del gem.data_mem[ds]
    ## end tiling

    ## set up image segments
    if setu['dsf_aot_estimate'] == 'segmented':
        segment_data = {}
        rhot_ds = [ds for ds in gem.datasets if 'rhot_' in ds]
        finite_mask = np.isfinite(gem.data(rhot_ds[0]))
        del rhot_ds
        segment_mask = skimage.measure.label(finite_mask)
        segments = np.unique(segment_mask)

        ## find and label segments
        for segment in segments:
            #if segment == 0: continue
            seg_sub = np.where((segment_mask == segment) & (finite_mask))
            #if len(seg_sub[0]) == 0: continue
            if len(seg_sub[0]) < max(1, setu['dsf_minimum_segment_size']):
                if setu['verbosity'] > 4: print('Skipping segment of {} pixels'.format(len(seg_sub[0])))
                continue
            segment_data[segment] = {'segment': segment, 'sub': seg_sub}

        if len(segment_data) <= 1:
            print('Image segmentation only found {} segments'.format(len(segment_data)))
            print('Proceeding with dsf_aot_estimate=fixed')
            setu['dsf_aot_estimate'] = 'fixed'
        else:
            if setu['verbosity'] > 3: print('Found {} segments'.format(len(segment_data)))
            for segment in segment_data:
                if setu['verbosity'] > 4: print('Segment {}/{}: {} pixels'.format(segment, len(segment_data), len(segment_data[segment]['sub'][0])))
            ## convert geometry and ancillary data
            for ds in geom_ds:
                if len(np.atleast_1d(gem.data(ds)))>1: ## if not fixed geometry
                    gem.data_mem['{}_segmented'.format(ds)] = [np.nanmean(gem.data(ds)[segment_data[segment]['sub']]) for segment in segment_data]
                else:
                    gem.data_mem['{}_segmented'.format(ds)] = [1.0 * gem.data(ds) for segment in segment_data]
                gem.data_mem['{}_segmented'.format(ds)] = np.asarray(gem.data_mem['{}_segmented'.format(ds)]).flatten()
    ## end segmenting

    if (not setu['resolved_geometry']) & (setu['dsf_aot_estimate'] != 'tiled'): use_revlut = False
    if setu['dsf_aot_estimate'] in ['fixed', 'segmented']: use_revlut = False
    if hyper: use_revlut = False

    ## set LUT dimension parameters to correct shape if resolved processing
    if (use_revlut) & (per_pixel_geometry) & (setu['dsf_aot_estimate'] == 'resolved'):
        for ds in geom_ds:
            if len(np.atleast_1d(gem.data(ds)))!=1: continue
            print('Reshaping {} to {}x{} pixels for resolved processing'.format(ds, gem.gatts['data_dimensions'][0], gem.gatts['data_dimensions'][1]))
            gem.data_mem[ds] = np.repeat(gem.data_mem[ds], gem.gatts['data_elements']).reshape(gem.gatts['data_dimensions'])


    ## set output directory
    if setu['output'] is not None:
        output_ = setu['output']
    else:
        output_ = os.path.dirname(gemf)

    ## write settings
    settings_file = '{}/acolite_run_{}_l2r_settings.txt'.format(output_,setu['runid'])
    ac.acolite.settings.write(settings_file, setu)

    ## do not allow LUT boundaries ()
    left, right = np.nan, np.nan
    if setu['dsf_allow_lut_boundaries']:
        left, right = None, None

    ## setup output file
    ofile = None
    if output_file:
        if target_file is None:
            ofile = gemf.replace('_L1R', '_L2R')
            #if ('output' in setu) & (output is None): output = setu['output']
            if output is None: output = output_
            ofile = '{}/{}'.format(output, os.path.basename(ofile))
        else:
            ofile = '{}'.format(target_file)

        gemo = ac.gem.gem(ofile, new=True)
        gemo.gatts = {k: gem.gatts[k] for k in gem.gatts}
        gemo.nc_projection = gem.nc_projection
        gemo.gatts['acolite_file_type'] = 'L2R'
        gemo.gatts['ofile'] = ofile

        gemo.bands = gem.bands
        gemo.verbosity = setu['verbosity']
        gemo.gatts['acolite_version'] = ac.version

        ## add settings to gatts
        for k in setu:
            if k in gem.gatts: continue
            if setu[k] in [True, False]:
                gemo.gatts[k] = str(setu[k])
            else:
                gemo.gatts[k] = setu[k]

        ## copy datasets from inputfile
        copy_rhot = False
        copy_datasets = []
        if setu['copy_datasets'] is not None: copy_datasets += setu['copy_datasets']
        if setu['output_bt']: copy_datasets += [ds for ds in gem.datasets if ds[0:2] == 'bt']
        if setu['output_xy']: copy_datasets += ['x', 'y']

        if len(copy_datasets) > 0:
            ## copy rhot all from L1R
            if 'rhot_*' in copy_datasets:
                copy_datasets.remove('rhot_*')
                copy_rhot = True
            ## copy datasets to L2R
            for ds in copy_datasets:
                if (ds not in gem.datasets):
                    if setu['verbosity'] > 2: print('{} not found in {}'.format(ds, gemf))
                    continue
                if setu['verbosity'] > 1: print('Writing {}'.format(ds))
                cdata, catts = gem.data(ds, attributes=True)
                gemo.write(ds, cdata, ds_att=catts)
                del cdata, catts

        ## write dem
        if setu['dem_pressure_write']:
            for k in ['dem', 'dem_pressure']:
                if k in gem.data_mem:
                    gemo.write(k, gem.data_mem[k])
                    gem.data_mem[k] = None

    t0 = time.time()
    print('Loading LUTs {}'.format(setu['luts']))
    ## load reverse lut romix -> aot
    if use_revlut: revl = ac.aerlut.reverse_lut(sensor_lut, par = ipar, \
                                                rsky_lut = setu['dsf_interface_lut'], base_luts=setu['luts'])
    ## load aot -> atmospheric parameters lut
    ## QV 2022-04-04 interface reflectance is always loaded since we include wind in the interpolation below
    ## not necessary for runs with par == romix, to be fixed
    lutdw = ac.aerlut.import_luts(add_rsky = True, par = ipar, sensor=None if hyper else sensor_lut,
                                  rsky_lut = setu['dsf_interface_lut'],
                                  base_luts = setu['luts'], pressures = setu['luts_pressures'],
                                  reduce_dimensions=setu['luts_reduce_dimensions'])
    luts = list(lutdw.keys())
    print('Loading LUTs took {:.1f} s'.format(time.time()-t0))

    ## #####################
    ## dark spectrum fitting
    if (ac_opt == 'dsf'):
        ## user supplied aot
        if (setu['dsf_fixed_aot'] is not None):
            aot_lut = None
            for li, lut in enumerate(luts):
                if lut == setu['dsf_fixed_lut']:
                    aot_lut = np.array(li)
                    aot_lut.shape+=(1,1) ## make 1,1 dimensions
            if aot_lut is None:
                print('LUT {} not recognised'.format(setu['dsf_fixed_lut']))
            aot_sel_lut = luts[aot_lut[0][0]]
            aot_sel_par = None

            ## ancillary aot can be 2D
            if setu['dsf_aot_estimate'].startswith('ancillary'):
                if setu['dsf_aot_estimate'] == 'ancillary_fixed':
                    aot_sel = np.array(np.nanmean(aer_aot))
                else:
                    aot_sel = aer_aot * 1.0
                del aer_aot
            else:
                aot_sel = np.array(float(setu['dsf_fixed_aot']))
                aot_sel.shape += (1,1)

            if len(np.atleast_1d(aot_sel)) > 1:
                print('User specified aot with shape {}x{} and model {}'.format(aot_sel.shape[0], aot_sel.shape[1], aot_sel_lut))
            else:
                aot_sel.shape += (1,1)
                print('User specified aot {:.3f} and model {}'.format(aot_sel[0][0], aot_sel_lut))

            ## set to fixed for processing
            setu['dsf_aot_estimate'] = 'fixed'
            ## geometry key '' if using resolved, otherwise '_mean' or '_tiled'
            gk = '' if use_revlut else '_mean'
        ## image derived aot
        else:
            if setu['dsf_spectrum_option'] not in ['darkest', 'percentile', 'intercept']:
                print('dsf_spectrum_option {} not configured, falling back to darkest'.format(setu['dsf_spectrum_option']))
                setu['dsf_spectrum_option'] = 'darkest'

            ## run through bands to get aot
            aot_bands = []
            aot_dict = {}
            dsf_rhod = {}
            dsf_location = {}
            for bi, b in enumerate(gem.bands):
                ## test if band can or should be used in DSF
                ## skip band if listed in dsf_exlude_bands
                if (str(b) in setu['dsf_exclude_bands']):
                    if setu['verbosity'] > 5: print('Skipping band {} ({}) as it is in dsf_exclude_bands: {}'.format(b, gem.bands[b]['rhot_ds'], setu['dsf_exclude_bands']))
                    continue
                ## skip bands according to configuration
                if (gem.bands[b]['wave_nm'] < setu['dsf_wave_range'][0]):
                    if setu['verbosity'] > 5: print('Skipping band {} ({}) as wavelength < dsf_wave_range[0]: {:.1f} < {}'.format(b, gem.bands[b]['rhot_ds'], gem.bands[b]['wave_nm'], setu['dsf_wave_range'][0]))
                    continue
                if (gem.bands[b]['wave_nm'] > setu['dsf_wave_range'][1]):
                    if setu['verbosity'] > 5: print('Skipping band {} ({}) as wavelength > dsf_wave_range[1]: {:.1f} > {}'.format(b, gem.bands[b]['rhot_ds'], gem.bands[b]['wave_nm'], setu['dsf_wave_range'][1]))
                    continue
                ## skip band for aot computation
                if gem.bands[b]['tt_gas'] < setu['min_tgas_aot']:
                    if setu['verbosity'] > 5: print('Skipping band {} ({}) as tt_gas < min_tgas_aot: {:.3f} < {:.3f}'.format(b, gem.bands[b]['rhot_ds'], gem.bands[b]['tt_gas'],setu['min_tgas_aot']))
                    continue
                ## skip band if either rhot_ds or tt_gas are missing from attributes
                if ('rhot_ds' not in gem.bands[b]) or ('tt_gas' not in gem.bands[b]):
                    if setu['verbosity'] > 5: print('Skipping band {} ({}) as rhot_ds or tt_gas is missing from attributes'.format(b, gem.bands[b]['rhot_ds']))
                    continue
                ## skip band if rhot data not in datasets
                if gem.bands[b]['rhot_ds'] not in gem.datasets:
                    if setu['verbosity'] > 5: print('Skipping band {} ({}) as {} not in datasets: {}'.format(b, gem.bands[b]['rhot_ds'], gem.bands[b]['rhot_ds'], gem.datasets))
                    continue
                ## end test if band can or should be used in DSF

                if setu['verbosity'] > 1: print('Running AOT estimation for band {} ({})'.format(b, gem.bands[b]['rhot_ds']))

                ## extract data
                band_data = gem.data(gem.bands[b]['rhot_ds'])*1.0
                band_shape = band_data.shape

                ## resample to reduce staircase effect
                ## apply filter for SNBC
                if (setu['sensor_noise_bias_correction']) & ('sensor_noise' in gem.bands[b]):
                    if 'sensor_noise_bias_correction_resampling' in setu:
                        if setu['sensor_noise_bias_correction_resampling'] is None:
                            n_noise = 1
                        else:
                            n_noise = int(setu['sensor_noise_bias_correction_resampling'])
                            if n_noise < 1:
                                n_noise = 1
                                print('Setting sensor_noise_bias_correction_resampling={}'.format(n_noise))
                            else:
                                print('Using sensor_noise_bias_correction_resampling={}'.format(n_noise))
                    else:
                        n_noise = 1

                    ## do resampling if n_noise > 1
                    if n_noise > 1:
                        d_ = band_data * 1.0
                        ## subset to number of pixels divisible by n_noise
                        m1 = d_.shape[0] - (d_.shape[0] % n_noise)
                        m2 = d_.shape[1] - (d_.shape[1] % n_noise)
                        d_ = d_[0:m1, 0:m2]

                        ## resampled size
                        n1 = int(d_.shape[0] / n_noise)
                        n2 = int(d_.shape[1] / n_noise)
                        d_2 = d_.reshape(n_noise, n_noise, n1, n2)

                        ## compute average over nxn pixels
                        d_mean = d_2.mean(axis = (0,1))
                        band_data = d_mean * 1.0
                ## end apply filter for SNBC

                ## compute mask
                valid = np.isfinite(band_data)*(band_data>0)
                mask = valid == False

                ## apply TOA filter
                if setu['dsf_filter_rhot']:
                    if setu['verbosity'] > 1: print('Filtered {} using {}th percentile in {}x{} pixel box'.format(gem.bands[b]['rhot_ds'],
                                                    setu['dsf_filter_percentile'], setu['dsf_filter_box'][0], setu['dsf_filter_box'][1]))
                    band_data[mask] = np.nanmedian(band_data) ## fill mask with median
                    #band_data = scipy.ndimage.median_filter(band_data, size=setu['dsf_filter_box'])
                    band_data = scipy.ndimage.percentile_filter(band_data, setu['dsf_filter_percentile'], size=setu['dsf_filter_box'])
                    band_data[mask] = np.nan
                band_sub = np.where(valid)
                del valid, mask

                ## geometry key '' if using resolved, otherwise '_mean' or '_tiled'
                gk = ''

                ## fixed path reflectance
                if setu['dsf_aot_estimate'] == 'fixed':
                    band_data_copy = band_data * 1.0
                    if setu['dsf_spectrum_option'] == 'darkest':
                        band_data = np.array((np.nanpercentile(band_data[band_sub], 0)))
                    if setu['dsf_spectrum_option'] == 'percentile':
                        band_data = np.array((np.nanpercentile(band_data[band_sub], setu['dsf_percentile'])))
                    if setu['dsf_spectrum_option'] == 'intercept':
                        band_data = ac.shared.intercept(band_data[band_sub], setu['dsf_intercept_pixels'])
                    band_data.shape+=(1,1) ## make 1,1 dimensions
                    gk='_mean'
                    dark_pixel_location = np.where(band_data_copy <= band_data)
                    dark_pixel_location_x, dark_pixel_location_y = -1, -1
                    if len(dark_pixel_location[0]) != 0: ## store brightest pixel lower than band data
                        pixid = np.argsort((band_data_copy - band_data)[dark_pixel_location])
                        dark_pixel_location_x = dark_pixel_location[1][pixid[-1]]
                        dark_pixel_location_y = dark_pixel_location[0][pixid[-1]]
                    dsf_location[b] = (dark_pixel_location_x, dark_pixel_location_y)
                    if setu['verbosity'] > 2: print(b, setu['dsf_spectrum_option'], '{:.3f}'.format(float(band_data[0,0])))

                ## tiled path reflectance
                elif setu['dsf_aot_estimate'] == 'tiled':
                    gk = '_tiled'

                    ## tile this band data
                    tile_data = np.zeros((tiles[-1][0]+1, tiles[-1][1]+1), dtype=np.float32) + np.nan
                    for t in range(len(tiles)):
                        ti, tj, subti, subtj = tiles[t]
                        tsub = band_data[subti[0]:subti[1], subtj[0]:subtj[1]]
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
                        del tsub

                    ## fill nan tiles with closest values
                    ind = scipy.ndimage.distance_transform_edt(np.isnan(tile_data), return_distances=False, return_indices=True)
                    band_data = tile_data[tuple(ind)]
                    del tile_data, ind

                ## image is segmented based on input vector mask
                elif setu['dsf_aot_estimate'] == 'segmented':
                    gk = '_segmented'
                    if setu['dsf_spectrum_option'] == 'darkest':
                        band_data = np.array([np.nanpercentile(band_data[segment_data[segment]['sub']], 0)[0] for segment in segment_data])
                    if setu['dsf_spectrum_option'] == 'percentile':
                        band_data = np.array([np.nanpercentile(band_data[segment_data[segment]['sub']], setu['dsf_percentile'])[0] for segment in segment_data])
                    if setu['dsf_spectrum_option'] == 'intercept':
                        band_data = np.array([ac.shared.intercept(band_data[segment_data[segment]['sub']], setu['dsf_intercept_pixels'])  for segment in segment_data])
                    band_data.shape+=(1,1) ## make 2 dimensions

                ## resolved per pixel dsf
                elif setu['dsf_aot_estimate'] == 'resolved':
                    if not setu['resolved_geometry']: gk = '_mean'
                else:
                    print('DSF option {} not configured'.format(setu['dsf_aot_estimate']))
                    continue
                del band_sub

                band_sub = np.where(np.isfinite(band_data))
                if len(band_sub[0]) == 0:
                    print('No valid TOA data for band {}'.format(b))
                    continue

                ## do noise bias correction
                if (setu['sensor_noise_bias_correction']) & ('sensor_noise' in gem.bands[b]):
                    if gem.bands[b]['sensor_noise'] > 0:
                        if setu['sensor_noise_bias_correction_sigma_factor'] is not None:
                            print('Performing sensor noise bias correction for band {} with {} x {}'.format(b, setu['sensor_noise_bias_correction_sigma_factor'], gem.bands[b]['sensor_noise']))
                            noise_offset = gem.bands[b]['sensor_noise'] * setu['sensor_noise_bias_correction_sigma_factor']
                        else:
                            if setu['dsf_spectrum_option'] == 'percentile':
                                zf = scipy.stats.norm.ppf(1-(setu['dsf_percentile'])/100)
                                print('Using sigma factor of {:.3f} based on dsf_percentile = {}'.format(zf, setu['dsf_percentile']))
                            else:
                                zf = 2
                                print('Warning: Use of sensor_noise_bias_correction without sigma factor is recommended for dsf_spectrum_option=percentile')
                                print('Warning: Using default sigma factor = {:.3f}'.format(zf))
                            print('Performing sensor noise bias correction for band {} with {:.3f} x {}'.format(b, zf, gem.bands[b]['sensor_noise']))
                            noise_offset = gem.bands[b]['sensor_noise'] * zf

                        ## reduce noise offset by 1/n
                        noise_offset *= 1/n_noise
                        print('Performing sensor noise bias correction with noise reduction of 1/{}'.format(n_noise))

                        if setu['sensor_noise_bias_correction_sun_zenith']:
                            szaf = np.cos(np.radians(gem.data_mem['sza'+gk][band_sub]))
                            print('Performing sensor noise bias correction with sun zenith angle factor applied (mean = {:.3f})'.format(1/np.nanmean(szaf)))
                            band_data += noise_offset / szaf
                            del szaf
                        else:
                            band_data += noise_offset
                ## end noise bias correction

                ## do gas correction
                band_data[band_sub] /= gem.bands[b]['tt_gas']

                ## store rhod
                if setu['dsf_aot_estimate'] in ['fixed', 'tiled', 'segmented']:
                    dsf_rhod[b] = band_data

                ## use band specific geometry if available
                gk_raa = '{}'.format(gk)
                gk_vza = '{}'.format(gk)
                if 'raa_{}'.format(gem.bands[b]['wave_name']) in gem.datasets:
                    gk_raa = '_{}'.format(gem.bands[b]['wave_name'])+gk_raa
                if 'vza_{}'.format(gem.bands[b]['wave_name']) in gem.datasets:
                    gk_vza = '_{}'.format(gem.bands[b]['wave_name'])+gk_vza

                ## compute aot
                aot_band = {}
                for li, lut in enumerate(luts):
                    rhot_aot = None ## set rhot_aot to None for next LUT
                    aot_band[lut] = np.zeros(band_data.shape, dtype=np.float32)+np.nan
                    t0 = time.time()

                    ## reverse lut interpolates rhot directly to aot
                    if use_revlut:
                        if len(revl[lut]['rgi'][b].grid) == 5:
                            aot_band[lut][band_sub] = revl[lut]['rgi'][b]((gem.data_mem['pressure'+gk][band_sub],
                                                                               gem.data_mem['raa'+gk_raa][band_sub],
                                                                               gem.data_mem['vza'+gk_vza][band_sub],
                                                                               gem.data_mem['sza'+gk][band_sub],
                                                                               band_data[band_sub]))
                        else:
                            aot_band[lut][band_sub] = revl[lut]['rgi'][b]((gem.data_mem['pressure'+gk][band_sub],
                                                                               gem.data_mem['raa'+gk_raa][band_sub],
                                                                               gem.data_mem['vza'+gk_vza][band_sub],
                                                                               gem.data_mem['sza'+gk][band_sub],
                                                                               gem.data_mem['wind'+gk][band_sub],
                                                                               band_data[band_sub]))
                        # mask out of range aot
                        aot_band[lut][aot_band[lut]<=revl[lut]['minaot']]=np.nan
                        aot_band[lut][aot_band[lut]>=revl[lut]['maxaot']]=np.nan

                        ## replace nans with closest aot
                        if (setu['dsf_aot_fillnan']): aot_band[lut] = ac.shared.fillnan(aot_band[lut])

                    ## standard lut interpolates rhot to results for different aot values
                    else:
                        ## get rho path for lut steps in aot
                        if hyper:
                            # get modeled rhot for each wavelength
                            if rhot_aot is None:
                                ## set up array to store modeled rhot
                                rhot_aot = np.zeros((len(lutdw[lut]['meta']['tau']), \
                                                     len(lutdw[lut]['meta']['wave']), \
                                                     len(gem.data_mem['pressure'+gk].flatten())))

                                ## compute rhot for range of aot
                                for ai, aot in enumerate(lutdw[lut]['meta']['tau']):
                                    for pi in range(rhot_aot.shape[2]):
                                        tmp = lutdw[lut]['rgi']((gem.data_mem['pressure'+gk].flatten()[pi],
                                                                  lutdw[lut]['ipd'][par],
                                                                  lutdw[lut]['meta']['wave'],
                                                                  gem.data_mem['raa'+gk_raa].flatten()[pi],
                                                                  gem.data_mem['vza'+gk_vza].flatten()[pi],
                                                                  gem.data_mem['sza'+gk].flatten()[pi],
                                                                  gem.data_mem['wind'+gk].flatten()[pi], aot))
                                        ## store current result
                                        rhot_aot[ai,:,pi] = tmp.flatten()
                                if setu['verbosity'] > 4: print('Shape of modeled rhot: {}'.format(rhot_aot.shape))
                            ## resample modeled results to current band
                            tmp = ac.shared.rsr_convolute_nd(rhot_aot, lutdw[lut]['meta']['wave'], rsrd['rsr'][b]['response'], rsrd['rsr'][b]['wave'], axis=1)

                            ## interpolate rho path to observation
                            aotret = np.zeros(aot_band[lut][band_sub].flatten().shape)
                            ## interpolate to observed rhot
                            for ri, crho in enumerate(band_data.flatten()):
                                aotret[ri] = np.interp(crho, tmp[:,ri], lutdw[lut]['meta']['tau'], left=left, right=right)
                            if setu['verbosity'] > 4: print('Shape of computed aot: {}'.format(aotret.shape))
                            if setu['verbosity'] > 5: print('Computed aot: {}'.format(aotret))
                            aot_band[lut][band_sub] = aotret.reshape(aot_band[lut][band_sub].shape)
                        else:
                            if len(gem.data_mem['pressure'+gk]) > 1:
                                for gki in range(len(gem.data_mem['pressure'+gk])):
                                    tmp = lutdw[lut]['rgi'][b]((gem.data_mem['pressure'+gk][gki],
                                                                lutdw[lut]['ipd'][par],
                                                                gem.data_mem['raa'+gk_raa][gki],
                                                                gem.data_mem['vza'+gk_vza][gki],
                                                                gem.data_mem['sza'+gk][gki],
                                                                gem.data_mem['wind'+gk][gki], lutdw[lut]['meta']['tau']))
                                    tmp = tmp.flatten()
                                    aot_band[lut][gki] = np.interp(band_data[gki], tmp, lutdw[lut]['meta']['tau'], left=left, right=right)
                            else:
                                tmp = lutdw[lut]['rgi'][b]((gem.data_mem['pressure'+gk],
                                                            lutdw[lut]['ipd'][par],
                                                            gem.data_mem['raa'+gk_raa],
                                                            gem.data_mem['vza'+gk_vza],
                                                            gem.data_mem['sza'+gk],
                                                            gem.data_mem['wind'+gk], lutdw[lut]['meta']['tau']))
                                tmp = tmp.flatten()

                                ## interpolate rho path to observation
                                aot_band[lut][band_sub] = np.interp(band_data[band_sub], tmp, lutdw[lut]['meta']['tau'], left=left, right=right)

                    ## mask minimum/maximum tile aots
                    if setu['dsf_aot_estimate'] == 'tiled':
                        aot_band[lut][aot_band[lut]<setu['dsf_min_tile_aot']]=np.nan
                        aot_band[lut][aot_band[lut]>setu['dsf_max_tile_aot']]=np.nan

                    tel = time.time()-t0
                    if setu['verbosity'] > 1: print('{}/B{} {} took {:.3f}s ({})'.format(sensor_lut, b, lut, tel, 'RevLUT' if use_revlut else 'StdLUT'))

                ## store current band results
                aot_dict[b] = aot_band
                aot_bands.append(b)
                del band_data, band_sub, aot_band

            ## test if valid data could be extracted
            if len(aot_bands) == 0:
                print('No valid data found for aot retrieval with current settings.')
                gemo.close()
                print('Deleting {}.'.format(ofile))
                os.remove(ofile)
                return()

            ## get min aot per pixel
            aot_stack = {}
            for li, lut in enumerate(luts):
                aot_band_list = []
                ## stack aot for this lut
                for bi, b in enumerate(aot_bands):
                    if b not in aot_dict: continue
                    aot_band_list.append(b)
                    if lut not in aot_stack:
                        aot_stack[lut] = {'all': np.atleast_3d(aot_dict[b][lut])}
                    else:
                        aot_stack[lut]['all'] = np.dstack((aot_stack[lut]['all'], aot_dict[b][lut]))
                aot_stack[lut]['band_list'] = aot_band_list

                ## sort aot per pixel
                tmp = np.argsort(aot_stack[lut]['all'], axis=2)
                ay, ax = np.meshgrid(np.arange(tmp.shape[1]), np.arange(tmp.shape[0]))

                ## identify number of bands
                if setu['dsf_nbands']<2: setu['dsf_nbands'] = 2
                if setu['dsf_nbands']>tmp.shape[2]: setu['dsf_nbands'] = tmp.shape[2]
                if setu['dsf_nbands_fit']<2: setu['dsf_nbands_fit'] = 2
                if setu['dsf_nbands_fit']>tmp.shape[2]: setu['dsf_nbands_fit'] = tmp.shape[2]

                ## get minimum or average aot
                if setu['dsf_aot_compute'] in ['mean', 'median']:
                    print('Using dsf_aot_compute = {}'.format(setu['dsf_aot_compute']))
                    ## stack n lowest bands
                    for ai in range(setu['dsf_nbands']):
                        if ai == 0:
                            tmp_aot = aot_stack[lut]['all'][ax, ay, tmp[ax,ay,ai]] * 1.0
                        else:
                            tmp_aot = np.dstack((tmp_aot, aot_stack[lut]['all'][ax, ay, tmp[ax,ay,ai]] * 1.0))
                    ## compute mean over stack
                    if setu['dsf_aot_compute'] == 'mean': aot_stack[lut]['aot'] = np.nanmean(tmp_aot, axis=2)
                    if setu['dsf_aot_compute'] == 'median': aot_stack[lut]['aot'] = np.nanmedian(tmp_aot, axis=2)
                    if setu['dsf_aot_estimate'] == 'fixed': print('Using dsf_aot_compute = {} {} aot = {:.3f}'.format(setu['dsf_aot_compute'], lut, float(aot_stack[lut]['aot'].flatten())))
                    tmp_aot = None
                else:
                    aot_stack[lut]['aot'] = aot_stack[lut]['all'][ax,ay,tmp[ax,ay,0]] #np.nanmin(aot_stack[lut]['all'], axis=2)

                ## if minimum for fixed retrieval is nan, set it to 0.01
                if setu['dsf_aot_estimate'] == 'fixed':
                    if np.isnan(aot_stack[lut]['aot']): aot_stack[lut]['aot'][0][0] = 0.01
                aot_stack[lut]['mask'] = ~np.isfinite(aot_stack[lut]['aot'])

                ## apply percentile filter
                if (setu['dsf_filter_aot']) & (setu['dsf_aot_estimate'] == 'resolved'):
                    aot_stack[lut]['aot'] = \
                        scipy.ndimage.percentile_filter(aot_stack[lut]['aot'],
                                                        setu['dsf_filter_percentile'],
                                                        size=setu['dsf_filter_box'])
                ## apply gaussian kernel smoothing
                if (setu['dsf_smooth_aot']) & (setu['dsf_aot_estimate'] == 'resolved'):
                    ## for gaussian smoothing of aot
                    aot_stack[lut]['aot'] = scipy.ndimage.gaussian_filter(aot_stack[lut]['aot'], setu['dsf_smooth_box'], order=0, mode='nearest')

                ## mask aot
                aot_stack[lut]['aot'][aot_stack[lut]['mask']] = np.nan

                ## store bands for fitting rmsd
                for bbi in range(setu['dsf_nbands_fit']):
                    aot_stack[lut]['b{}'.format(bbi+1)] = tmp[:,:,bbi].astype(int)#.astype(float)
                    aot_stack[lut]['b{}'.format(bbi+1)][aot_stack[lut]['mask']] = -1

                if setu['dsf_model_selection'] == 'min_dtau':
                    ## array idices
                    aid = np.indices(aot_stack[lut]['all'].shape[0:2])
                    ## abs difference between first and second band tau
                    aot_stack[lut]['dtau'] = np.abs(aot_stack[lut]['all'][aid[0,:],aid[1,:],tmp[:,:,0]]-\
                                                    aot_stack[lut]['all'][aid[0,:],aid[1,:],tmp[:,:,1]])
                ## remove sorted indices
                tmp = None

            ## use only model in luts
            if len(luts) == 1:
                lut = luts[0]
                fit_bands = ['b{}'.format(bbi+1) for bbi in range(setu['dsf_nbands_fit'])]
                if setu['verbosity'] > 1: print('Using only model: {}'.format(lut))

                aot_lut = np.zeros(aot_stack[lut]['aot'].shape, dtype=np.float32).astype(int)
                aot_lut[aot_stack[lut]['mask']] = -1
                aot_sel = aot_stack[lut]['aot'] * 1.0
                aot_sel_par = aot_stack[lut]['aot'] * np.nan
                aot_sel_lut = '{}'.format(lut)
                if setu['dsf_aot_estimate'] == 'fixed':
                    aot_sel_bands = [aot_stack[lut]['{}'.format(bb)][0][0] for bb in fit_bands]
            ## select model based on min rmsd for 2 bands
            else:
                if setu['verbosity'] > 1: print('Choosing best fitting model: {} ({} bands)'.format(setu['dsf_model_selection'], setu['dsf_nbands']))

                ## run through model results, get rhod and rhop for n lowest bands
                for li, lut in enumerate(luts):
                    ## select model based on minimum rmsd between n best fitting bands
                    if setu['dsf_model_selection'] == 'min_drmsd':
                        if setu['verbosity'] > 1: print('Computing RMSD for model {}'.format(lut))
                        rhop_f = np.zeros((aot_stack[lut]['b1'].shape[0],aot_stack[lut]['b1'].shape[1],setu['dsf_nbands_fit']), dtype=np.float32) + np.nan
                        rhod_f = np.zeros((aot_stack[lut]['b1'].shape[0],aot_stack[lut]['b1'].shape[1],setu['dsf_nbands_fit']), dtype=np.float32) + np.nan
                        for bi, b in enumerate(aot_bands):

                            ## use band specific geometry if available
                            gk_raa = '{}'.format(gk)
                            gk_vza = '{}'.format(gk)
                            if 'raa_{}'.format(gem.bands[b]['wave_name']) in gem.datasets:
                                gk_raa = '_{}'.format(gem.bands[b]['wave_name'])+gk_raa
                            if 'vza_{}'.format(gem.bands[b]['wave_name']) in gem.datasets:
                                gk_vza = '_{}'.format(gem.bands[b]['wave_name'])+gk_vza

                            ## run through two best fitting bands
                            fit_bands = ['b{}'.format(bbi+1) for bbi in range(setu['dsf_nbands_fit'])]
                            for ai, ab in enumerate(fit_bands):
                                aot_sub = np.where(aot_stack[lut][ab]==bi)
                                ## get rhod for current band
                                if (setu['dsf_aot_estimate'] == 'resolved'):
                                    rhod_f[aot_sub[0], aot_sub[1], ai] = gem.data(gem.bands[b]['rhot_ds'])[aot_sub]
                                elif (setu['dsf_aot_estimate'] == 'segmented'):
                                    rhod_f[aot_sub[0], aot_sub[1], ai] = dsf_rhod[b][aot_sub].flatten()
                                else:
                                    rhod_f[aot_sub[0], aot_sub[1], ai] = dsf_rhod[b][aot_sub]
                                ## get rho path for current band
                                if len(aot_sub[0]) > 0:
                                    if (use_revlut):
                                        xi = [gem.data_mem['pressure'+gk][aot_sub],
                                                          gem.data_mem['raa'+gk_raa][aot_sub],
                                                          gem.data_mem['vza'+gk_vza][aot_sub],
                                                          gem.data_mem['sza'+gk][aot_sub],
                                                          gem.data_mem['wind'+gk][aot_sub]]
                                    else:
                                        xi = [gem.data_mem['pressure'+gk],
                                                          gem.data_mem['raa'+gk_raa],
                                                          gem.data_mem['vza'+gk_vza],
                                                          gem.data_mem['sza'+gk],
                                                          gem.data_mem['wind'+gk]]
                                    if hyper:
                                        ## get hyperspectral results and resample to band
                                        if len(aot_stack[lut]['aot'][aot_sub]) == 1:
                                            if len(xi[0]) == 0:
                                                res_hyp = lutdw[lut]['rgi']((xi[0], lutdw[lut]['ipd'][par], lutdw[lut]['meta']['wave'],
                                                                                        xi[1], xi[2], xi[3], xi[4], aot_stack[lut]['aot'][aot_sub]))
                                            else: ## if more resolved geometry
                                                res_hyp = lutdw[lut]['rgi']((xi[0][aot_sub], lutdw[lut]['ipd'][par], lutdw[lut]['meta']['wave'],
                                                                             xi[1][aot_sub], xi[2][aot_sub], xi[3][aot_sub], xi[4][aot_sub], aot_stack[lut]['aot'][aot_sub]))
                                            rhop_f[aot_sub[0], aot_sub[1], ai] = ac.shared.rsr_convolute_nd(res_hyp.flatten(), lutdw[lut]['meta']['wave'], rsrd['rsr'][b]['response'], rsrd['rsr'][b]['wave'], axis=0)
                                        else:
                                            for iii in range(len(aot_stack[lut]['aot'][aot_sub])):
                                                if len(xi[0]) == 0:
                                                    res_hyp = lutdw[lut]['rgi']((xi[0], lutdw[lut]['ipd'][par], lutdw[lut]['meta']['wave'],
                                                                                            xi[1], xi[2], xi[3], xi[4], aot_stack[lut]['aot'][aot_sub][iii]))

                                                else: ## if more resolved geometry
                                                    res_hyp = lutdw[lut]['rgi']((xi[0].flatten()[iii], lutdw[lut]['ipd'][par], lutdw[lut]['meta']['wave'],
                                                                                            xi[1].flatten()[iii], xi[2].flatten()[iii], xi[3].flatten()[iii], xi[4].flatten()[iii], aot_stack[lut]['aot'][aot_sub][iii]))
                                                rhop_f[aot_sub[0][iii], aot_sub[1][iii], ai] = ac.shared.rsr_convolute_nd(res_hyp.flatten(), lutdw[lut]['meta']['wave'], rsrd['rsr'][b]['response'], rsrd['rsr'][b]['wave'], axis=0)
                                    else:
                                        if setu['dsf_aot_estimate'] == 'segmented':
                                            for gki in range(len(aot_sub[0])):
                                                rhop_f[aot_sub[0][gki], aot_sub[1][gki], ai] = lutdw[lut]['rgi'][b]((xi[0][aot_sub[0][gki]], lutdw[lut]['ipd'][par],
                                                                xi[1][aot_sub[0][gki]], xi[2][aot_sub[0][gki]],
                                                                xi[3][aot_sub[0][gki]], xi[4][aot_sub[0][gki]], aot_stack[lut]['aot'][aot_sub][gki]))

                                        else:
                                            rhop_f[aot_sub[0], aot_sub[1], ai] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd'][par],
                                                                                    xi[1], xi[2], xi[3], xi[4], aot_stack[lut]['aot'][aot_sub]))
                        ## rmsd for current bands
                        cur_sel_par = np.sqrt(np.nanmean(np.square((rhod_f-rhop_f)), axis=2))
                        if (setu['dsf_aot_estimate'] == 'fixed') & (setu['verbosity'] > 1): print('Computing RMSD for model {}: {:.4e}'.format(lut, cur_sel_par[0][0]))
                    ## end select with min RMSD

                    ## select model based on minimum delta tau between two lowest aot bands
                    if setu['dsf_model_selection'] == 'min_dtau':
                        cur_sel_par = aot_stack[lut]['dtau']
                    ## end select with min delta tau

                    ## store minimum info
                    if li == 0:
                        aot_lut = np.zeros(aot_stack[lut]['aot'].shape, dtype=np.float32).astype(int)
                        aot_lut[aot_stack[lut]['mask']] = -1
                        aot_sel = aot_stack[lut]['aot'] * 1.0
                        aot_sel_par = cur_sel_par * 1.0
                        if setu['dsf_aot_estimate'] == 'fixed':
                            aot_sel_lut = '{}'.format(lut)
                            aot_sel_bands = [aot_stack[lut]['{}'.format(bb)][0][0] for bb in fit_bands]
                    else:
                        aot_sub = np.where(cur_sel_par<aot_sel_par)
                        if len(aot_sub[0]) == 0: continue
                        aot_lut[aot_sub] = li
                        aot_sel[aot_sub] = aot_stack[lut]['aot'][aot_sub]*1.0
                        aot_sel_par[aot_sub] = cur_sel_par[aot_sub] * 1.0
                        if setu['dsf_aot_estimate'] == 'fixed':
                            aot_sel_lut = '{}'.format(lut)
                            aot_sel_bands = [aot_stack[lut]['{}'.format(bb)][0][0] for bb in fit_bands]

                rhod_f = None
                rhod_p = None

        ## check variable aot, use most common LUT
        if (setu['dsf_aot_estimate'] != 'fixed') & (setu['dsf_aot_most_common_model']) & (len(luts) > 1):
            print('Selecting most common model for processing.')
            n_aot = len(np.where(aot_lut != -1)[0])
            n_sel = 0
            for li, lut in enumerate(luts):
                sub = np.where(aot_lut == li)
                n_cur = len(sub[0])
                if n_cur == 0:
                    print('{}: {:.1f}%'.format(lut, 0))
                else:
                    print('{}: {:.1f}%: mean aot of subset = {:.2f}'.format(lut, 100*n_cur/n_aot, np.nanmean(aot_sel[sub])))
                if n_cur >= n_sel:
                    n_sel = n_cur
                    li_sel = li
                    aot_sel_lut = '{}'.format(lut)
            ## set selected model
            aot_lut[:] = li_sel
            aot_sel[:] = aot_stack[aot_sel_lut]['aot'][:]*1.0
            aot_sel_par[:] = np.nan # to do
            print('Selected {}, mean aot = {:.2f}'.format(aot_sel_lut, np.nanmean(aot_sel)))

        if (setu['dsf_aot_estimate'] == 'fixed') & (setu['verbosity'] > 1) & (aot_sel_par is not None):
            print('Selected model {}: aot {:.3f}, RMSD {:.2e}'.format(aot_sel_lut, aot_sel[0][0], aot_sel_par[0][0]))

    ### end dark_spectrum_fitting

    ## exponential
    elif ac_opt == 'exp':
        ## find bands to use
        exp_b1 = None
        exp_b1_diff = 1000
        exp_b2 = None
        exp_b2_diff = 1000
        exp_mask = None
        exp_mask_diff = 1000
        for b in gem.bands:
            sd = np.abs(gem.bands[b]['wave_nm'] - setu['exp_wave1'])
            if (sd < 100) & (sd < exp_b1_diff):
                exp_b1_diff = sd
                exp_b1 = b
                short_wv = gem.bands[b]['wave_nm']
            sd = np.abs(gem.bands[b]['wave_nm'] - setu['exp_wave2'])
            if (sd < 100) & (sd < exp_b2_diff):
                exp_b2_diff = sd
                exp_b2 = b
                long_wv = gem.bands[b]['wave_nm']
            sd = np.abs(gem.bands[b]['wave_nm'] - setu['exp_mask_wave'])
            if (sd < 100) & (sd < exp_mask_diff):
                exp_mask_diff = sd
                exp_mask = b
                mask_wv = gem.bands[b]['wave_nm']

        if (exp_b1 is None) or (exp_b2 is None):
            print('No suitable bands selected for EXP processing: exp_b1={}, exp_b2={}'.format(exp_b1, exp_b2))
            print('Please change wavelengths for exp_wave1={} and exp_wave2={} in your settings.'.format(setu['exp_wave1']. setu['exp_wave2']))
            return

        print('Selected bands {} and {} for EXP processing'.format(exp_b1, exp_b2))
        if (gem.bands[exp_b1]['rhot_ds'] not in gem.datasets) | (gem.bands[exp_b2]['rhot_ds'] not in gem.datasets):
            print('Selected bands are not available in {}'.format(gemf))
            if (gem.bands[exp_b1]['rhot_ds'] not in gem.datasets): print('EXP B1: {}'.format(gem.bands[exp_b1]['rhot_ds']))
            if (gem.bands[exp_b2]['rhot_ds'] not in gem.datasets): print('EXP B2: {}'.format(gem.bands[exp_b2]['rhot_ds']))
            return()

        ## determine processing option
        if (short_wv < 900) & (long_wv < 900):
            exp_option = 'red/NIR'
        elif (short_wv < 900) & (long_wv > 1500):
            exp_option = 'NIR/SWIR'
        else:
            exp_option = 'SWIR'

        ## read data
        exp_d1 = gem.data(gem.bands[exp_b1]['rhot_ds'])*1.0
        exp_d2 = gem.data(gem.bands[exp_b2]['rhot_ds'])*1.0

        if setu['resolved_geometry']:
            ## resolved
            xi = [gem.data_mem['pressure'],
                              gem.data_mem['raa'],
                              gem.data_mem['vza'],
                              gem.data_mem['sza'],
                              gem.data_mem['wind']]
        else:
            ## use mean geometry
            xi = [gem.data_mem['pressure'+'_mean'][0][0],
                  gem.data_mem['raa'+'_mean'][0][0],
                  gem.data_mem['vza'+'_mean'][0][0],
                  gem.data_mem['sza'+'_mean'][0][0],
                  gem.data_mem['wind'+'_mean'][0][0]]


        exp_lut = luts[0]
        exp_cwlim = 0.005
        exp_initial_epsilon = 1.0

        ## Rayleigh reflectance
        rorayl_b1 = lutdw[exp_lut]['rgi'][exp_b1]((xi[0], lutdw[exp_lut]['ipd'][par], xi[1], xi[2], xi[3], xi[4], 0.001))
        rorayl_b2 = lutdw[exp_lut]['rgi'][exp_b2]((xi[0], lutdw[exp_lut]['ipd'][par], xi[1], xi[2], xi[3], xi[4], 0.001))

        ## subtract Rayleigh reflectance
        exp_d1 -= rorayl_b1
        exp_d2 -= rorayl_b2

        ## compute mask
        if exp_mask == exp_b1:
            mask = exp_d1 >= setu['exp_mask_threshold']
        elif exp_mask == exp_b2:
            mask = exp_d2 >= setu['exp_mask_threshold']
        else:
            exp_dm = gem.data(gem.bands[exp_mask]['rhot_ds'])*1.0
            rorayl_mask = lutdw[exp_lut]['rgi'][exp_mask]((xi[0], lutdw[exp_lut]['ipd'][par], xi[1], xi[2], xi[3], xi[4], 0.001))
            exp_dm -= rorayl_mask
            mask = exp_dm >= setu['exp_mask_threshold']
            exp_dm = None
        print('Rayleigh corrected band {} used with threshold {:.3f} for clear water masking'.format(exp_mask, setu['exp_mask_threshold']))

        ## compute aerosol epsilon band ratio
        epsilon = exp_d1/exp_d2
        epsilon[np.where(mask)] = np.nan

        ## red/NIR option
        exp_fixed_epsilon = False
        if setu['exp_fixed_epsilon']: exp_fixed_epsilon = True

        if exp_option == 'red/NIR':
            print('Using similarity spectrum for red/NIR EXP')
            exp_fixed_epsilon = True

            ## Rayleigh transmittances in both bands
            dtotr_b1 = lutdw[exp_lut]['rgi'][exp_b1]((xi[0], lutdw[exp_lut]['ipd']['dtott'], xi[1], xi[2], xi[3], xi[4], 0.001))
            utotr_b1 = lutdw[exp_lut]['rgi'][exp_b1]((xi[0], lutdw[exp_lut]['ipd']['utott'], xi[1], xi[2], xi[3], xi[4], 0.001))
            dtotr_b2 = lutdw[exp_lut]['rgi'][exp_b2]((xi[0], lutdw[exp_lut]['ipd']['dtott'], xi[1], xi[2], xi[3], xi[4], 0.001))
            utotr_b2 = lutdw[exp_lut]['rgi'][exp_b2]((xi[0], lutdw[exp_lut]['ipd']['utott'], xi[1], xi[2], xi[3], xi[4], 0.001))
            tr_b1 = (dtotr_b1 * utotr_b1 * gem.bands[exp_b1]['tt_gas'])
            tr_b2 = (dtotr_b2 * utotr_b2 * gem.bands[exp_b2]['tt_gas'])
            ## get gamma
            if setu['exp_gamma'] is None:
                exp_gamma = tr_b1 / tr_b2
                print('Mean gamma: {:.2f}'.format(np.nanmean(exp_gamma)))
            else:
                exp_gamma = float(setu['exp_gamma'])
                print('User set gamma: {:.2f}'.format(exp_gamma))

            ## get alpha
            if setu['exp_alpha'] is None:
                ## import simspec
                simspec = ac.shared.similarity_read()
                ## convolute to sensor_o
                if setu['exp_alpha_weighted']:
                    ssd = ac.shared.rsr_convolute_dict(simspec['wave'], simspec['ave'], rsrd['rsr'])
                    exp_alpha = ssd[exp_b1]/ssd[exp_b2]
                ## or use closest bands
                else:
                    ssi0, ssw0 = ac.shared.closest_idx(simspec['wave'], gem.bands[exp_b1]['wave_mu'])
                    ssi1, ssw1 = ac.shared.closest_idx(simspec['wave'], gem.bands[exp_b2]['wave_mu'])
                    exp_alpha = simspec['ave'][ssi0]/simspec['ave'][ssi1]
                print('SimSpec alpha: {:.2f}'.format(exp_alpha))
            else:
                exp_alpha = float(setu['exp_alpha'])
                print('User set alpha: {:.2f}'.format(exp_alpha))

            ## first estimate of rhow to find clear waters
            exp_c1 = (exp_alpha/tr_b2)/(exp_alpha*exp_gamma-exp_initial_epsilon)
            exp_c2 = exp_initial_epsilon * exp_c1
            rhow_short = (exp_c1 * exp_d1) - (exp_c2 * exp_d2)

            ## additional masking for epsilon
            epsilon[(rhow_short < 0.) & (rhow_short > exp_cwlim)] = np.nan
            rhow_short = None
        elif exp_option == 'NIR/SWIR':
            print('Using NIR/SWIR EXP')
            exp_fixed_epsilon = True
            ## additional masking for epsilon
            mask2 = (exp_d2 < ((exp_d1+0.005) * 1.5) ) &\
                    (exp_d2 > ((exp_d1-0.005) * 0.8) ) &\
                    ((exp_d2 + 0.005)/exp_d1 > 0.8)
            epsilon[mask2] = np.nan
            mask2 = None
        elif exp_option == 'SWIR':
            print('Using SWIR EXP')
            if setu['exp_fixed_aerosol_reflectance']: exp_fixed_epsilon = True

        ## compute fixed epsilon
        if exp_fixed_epsilon:
            if setu['exp_epsilon'] is not None:
                epsilon = float(setu['exp_epsilon'])
            else:
                epsilon = np.nanpercentile(epsilon,setu['exp_fixed_epsilon_percentile'])

        ## determination of rhoam in long wavelength
        if exp_option == 'red/NIR':
            rhoam = (exp_alpha * exp_gamma * exp_d2 - exp_d1) / (exp_alpha * exp_gamma - epsilon)
        else:
            rhoam = exp_d2*1.0

        ## clear memory
        exp_d1,exp_d2 = None, None

        ## fixed rhoam?
        exp_fixed_rhoam = setu['exp_fixed_aerosol_reflectance']
        if exp_fixed_rhoam:
            rhoam = np.nanpercentile(rhoam,setu['exp_fixed_aerosol_reflectance_percentile'])
            print('{:.0f}th percentile rhoam ({:.0f} nm): {:.5f}'.format(setu['exp_fixed_aerosol_reflectance_percentile'], long_wv, rhoam))

        print('EXP band 1', setu['exp_wave1'], exp_b1, gem.bands[exp_b1]['rhot_ds'])
        print('EXP band 2', setu['exp_wave2'], exp_b2, gem.bands[exp_b2]['rhot_ds'])
        if exp_fixed_epsilon: print('Epsilon: {:.2f}'.format(epsilon))

        ## output data
        if setu['exp_output_intermediate']:
            if not exp_fixed_epsilon:   gemo.write('epsilon', epsilon)
            if not exp_fixed_rhoam: gemo.write('rhoam', rhoam)
    ## end exponential

    ## set up interpolator for tiled processing
    if (ac_opt == 'dsf') & (setu['dsf_aot_estimate'] == 'tiled'):
        xnew = np.linspace(0, tiles[-1][1], gem.gatts['data_dimensions'][1], dtype=np.float32)
        ynew = np.linspace(0, tiles[-1][0], gem.gatts['data_dimensions'][0], dtype=np.float32)

    ## store fixed aot in gatts
    if (ac_opt == 'dsf') & (setu['dsf_aot_estimate'] == 'fixed'):
        if aot_sel.shape == (1,1):
            gemo.gatts['ac_aot_550'] = aot_sel[0][0]
            gemo.gatts['ac_model'] = luts[aot_lut[0][0]]

        if setu['dsf_fixed_aot'] is None:
            ## store fitting parameter
            gemo.gatts['ac_fit'] = aot_sel_par[0][0]
            ## store bands used for DSF
            gemo.gatts['ac_bands'] = ','.join([str(b) for b in aot_stack[gemo.gatts['ac_model']]['band_list']])
            gemo.gatts['ac_nbands_fit'] = setu['dsf_nbands']
            ## store dark location
            for band in dsf_location:
                gemo.gatts['{}_xy_location'.format(band)] = dsf_location[band]
            ## store dark reflectance
            for band in dsf_rhod:
                gemo.gatts['{}_rhod'.format(band)] = dsf_rhod[band].flatten()

            for bbi, bn in enumerate(aot_sel_bands):
                gemo.gatts['ac_band{}_idx'.format(bbi+1)] = aot_sel_bands[bbi]
                gemo.gatts['ac_band{}'.format(bbi+1)] = aot_stack[gemo.gatts['ac_model']]['band_list'][aot_sel_bands[bbi]]

    ## store most common lut
    if (ac_opt == 'dsf') & (setu['dsf_aot_estimate'] != 'fixed') & (setu['dsf_aot_most_common_model']):
        gemo.gatts['ac_model'] = aot_sel_lut
        gemo.gatts['ac_aot_550_mean'] = np.nanmean(aot_sel)

    ## write aot to outputfile
    if (output_file) & (ac_opt == 'dsf') & (setu['dsf_write_aot_550']):
        ## reformat & save aot
        if setu['dsf_aot_estimate'] == 'fixed':
            if aot_sel.shape == (1,1):
                aot_out = np.repeat(aot_sel[0][0], gem.gatts['data_elements']).reshape(gem.gatts['data_dimensions'])
            else:
                aot_out = aot_sel * 1.0
        elif setu['dsf_aot_estimate'] == 'segmented':
            aot_out = np.zeros(gem.gatts['data_dimensions']) + np.nan
            for sidx, segment in enumerate(segment_data):
                aot_out[segment_data[segment]['sub']] = aot_sel[sidx]
        elif setu['dsf_aot_estimate'] == 'tiled':
            aot_out = ac.shared.tiles_interp(aot_sel, xnew, ynew, target_mask=None, smooth=setu['dsf_tile_smoothing'], kern_size=setu['dsf_tile_smoothing_kernel_size'], method=setu['dsf_tile_interp_method'])
        else:
            aot_out = aot_sel * 1.0
        ## write aot
        gemo.write('aot_550', aot_out)
        aot_out = None

    ## store ttot for glint correction
    ttot_all = {}

    ## allow use of per pixel geometry for fixed dsf
    if (per_pixel_geometry) & (setu['dsf_aot_estimate'] == 'fixed') & (setu['resolved_geometry']):
        use_revlut = True

    ## for ease of subsetting later, repeat single element datasets to the tile shape
    if (use_revlut) & (ac_opt == 'dsf') & (setu['dsf_aot_estimate'] != 'tiled'):
        for ds in geom_ds:
            if len(np.atleast_1d(gem.data(ds)))!=1: continue
            if setu['verbosity'] > 2: print('Reshaping {} to {}x{}'.format(ds, gem.gatts['data_dimensions'][0], gem.gatts['data_dimensions'][1]))
            gem.data_mem[ds] = np.repeat(gem.data_mem[ds], gem.gatts['data_elements']).reshape(gem.gatts['data_dimensions'])

    ## figure out cirrus bands
    if setu['cirrus_correction']:
        rho_cirrus = None

        ## use mean geometry to compute cirrus band Rayleigh
        xi = [gem.data_mem['pressure'+'_mean'][0][0],
              gem.data_mem['raa'+'_mean'][0][0],
              gem.data_mem['vza'+'_mean'][0][0],
              gem.data_mem['sza'+'_mean'][0][0],
              gem.data_mem['wind'+'_mean'][0][0]]

        ## compute Rayleigh reflectance for hyperspectral sensors
        if hyper:
            rorayl_hyp = lutdw[luts[0]]['rgi']((xi[0], lutdw[luts[0]]['ipd'][par],
                         lutdw[luts[0]]['meta']['wave'], xi[1], xi[2], xi[3], xi[4], 0.001)).flatten()

        ## find cirrus bands
        for bi, b in enumerate(gem.bands):
            if ('rhot_ds' not in gem.bands[b]): continue
            if gem.bands[b]['rhot_ds'] not in gem.datasets: continue
            if (gem.bands[b]['wave_nm'] < setu['cirrus_range'][0]): continue
            if (gem.bands[b]['wave_nm'] > setu['cirrus_range'][1]): continue

            ## compute Rayleigh reflectance
            if hyper:
                rorayl_cur = ac.shared.rsr_convolute_nd(rorayl_hyp, lutdw[luts[0]]['meta']['wave'], rsrd['rsr'][b]['response'], rsrd['rsr'][b]['wave'], axis=0)
            else:
                rorayl_cur = lutdw[luts[0]]['rgi'][b]((xi[0], lutdw[luts[0]]['ipd'][par], xi[1], xi[2], xi[3], xi[4], 0.001))

            ## cirrus reflectance = rho_t - rho_Rayleigh
            cur_data = gem.data(gem.bands[b]['rhot_ds']) - rorayl_cur

            if rho_cirrus is None:
                rho_cirrus = cur_data * 1.0
            else:
                rho_cirrus = np.dstack((rho_cirrus, cur_data))
            cur_data = None

        if rho_cirrus is None:
            setu['cirrus_correction'] = False
        else:
            ## compute mean from several bands
            if len(rho_cirrus.shape) == 3:
                rho_cirrus = np.nanmean(rho_cirrus, axis=2)
            ## write cirrus mean
            gemo.write('rho_cirrus', rho_cirrus)
    print('use_revlut', use_revlut)

    hyper_res = None
    ## compute surface reflectances
    for bi, b in enumerate(gem.bands):
        if ('rhot_ds' not in gem.bands[b]) or ('tt_gas' not in gem.bands[b]):
            if setu['verbosity'] > 2: print('Band {} at {} nm not in bands dataset'.format(b, gem.bands[b]['wave_name']))
            continue
        if gem.bands[b]['rhot_ds'] not in gem.datasets:
            if setu['verbosity'] > 2: print('Band {} at {} nm not in available rhot datasets'.format(b, gem.bands[b]['wave_name']))
            continue ## skip if we don't have rhot for a band that is in the RSR file

        ## temporary fix
        if (gem.bands[b]['wave_mu'] < 0.345) & (lutdw[luts[0]]['meta']['wave'][0] >= 0.34):
            if setu['verbosity'] > 2: print('Band {} at {} nm wavelength < 345 nm'.format(b, gem.bands[b]['wave_name']))
            continue ## skip if below LUT range

        dsi = gem.bands[b]['rhot_ds']
        dso = gem.bands[b]['rhos_ds']
        cur_data, cur_att = gem.data(dsi, attributes=True)

        ## store rhot in output file
        if copy_rhot:
            gemo.write(dsi, cur_data, ds_att = cur_att)

        if gem.bands[b]['tt_gas'] < setu['min_tgas_rho']:
            if setu['verbosity'] > 2: print('Band {} at {} nm has tgas < min_tgas_rho ({:.2f} < {:.2f})'.format(b, gem.bands[b]['wave_name'], gem.bands[b]['tt_gas'], setu['min_tgas_rho']))
            continue

        ## apply cirrus correction
        if setu['cirrus_correction']:
            g = setu['cirrus_g_vnir'] * 1.0
            if gem.bands[b]['wave_nm'] > 1000: g = setu['cirrus_g_swir'] * 1.0
            cur_data -= (rho_cirrus * g)

        t0 = time.time()
        if setu['verbosity'] > 1: print('Computing surface reflectance', b, gem.bands[b]['wave_name'], '{:.3f}'.format(gem.bands[b]['tt_gas']))

        ds_att = gem.bands[b]
        ds_att['wavelength']=ds_att['wave_nm']

        ## dark spectrum fitting
        if (ac_opt == 'dsf'):
            if setu['slicing']: valid_mask = np.isfinite(cur_data)

            ## shape of atmospheric datasets
            atm_shape = aot_sel.shape

            ## if path reflectance is resolved, but resolved geometry available
            if (use_revlut) & (setu['dsf_aot_estimate'] == 'fixed'):
                atm_shape = cur_data.shape
                gk = ''

            ## use band specific geometry if available
            gk_raa = '{}'.format(gk)
            gk_vza = '{}'.format(gk)
            if 'raa_{}'.format(gem.bands[b]['wave_name']) in gem.datasets:
                gk_raa = '_{}'.format(gem.bands[b]['wave_name'])+gk_raa
            if 'vza_{}'.format(gem.bands[b]['wave_name']) in gem.datasets:
                gk_vza = '_{}'.format(gem.bands[b]['wave_name'])+gk_vza

            romix = np.zeros(atm_shape, dtype=np.float32)+np.nan
            astot = np.zeros(atm_shape, dtype=np.float32)+np.nan
            dutott = np.zeros(atm_shape, dtype=np.float32)+np.nan
            if (setu['dsf_residual_glint_correction']) & (setu['dsf_residual_glint_correction_method']=='default'):
                ttot_all[b] = np.zeros(atm_shape, dtype=np.float32)+np.nan
            if (setu['output_ed']): dtott = np.zeros(atm_shape, dtype=np.float32)+np.nan

            for li, lut in enumerate(luts):
                ls = np.where(aot_lut == li)
                if len(ls[0]) == 0: continue
                ai = aot_sel[ls]

                ## resolved geometry with fixed path reflectance
                if (use_revlut) & (setu['dsf_aot_estimate'] == 'fixed'):
                    ls = np.where(cur_data)

                if (use_revlut):
                    xi = [gem.data_mem['pressure'+gk][ls],
                          gem.data_mem['raa'+gk_raa][ls],
                          gem.data_mem['vza'+gk_vza][ls],
                          gem.data_mem['sza'+gk][ls],
                          gem.data_mem['wind'+gk][ls]]
                else:
                    xi = [gem.data_mem['pressure'+gk],
                          gem.data_mem['raa'+gk_raa],
                          gem.data_mem['vza'+gk_vza],
                          gem.data_mem['sza'+gk],
                          gem.data_mem['wind'+gk]]
                    # subset to number of estimates made for this LUT
                    ## QV 2022-07-28 maybe not needed any more?
                    #if len(xi[0]) > 1:
                    #    xi = [[x[l] for l in ls[0]] for x in xi]

                if hyper:
                    ## compute hyper results and resample later
                    if hyper_res is None:
                        hyper_res = {}
                        for prm in [par, 'astot', 'dutott', 'ttot', 'dtott']:
                            if len(ai) == 1: ## fixed DSF
                                hyper_res[prm] = lutdw[lut]['rgi']((xi[0], lutdw[lut]['ipd'][prm],
                                                 lutdw[lut]['meta']['wave'], xi[1], xi[2], xi[3], xi[4], ai)).flatten()
                            else: ## tiled/resolved DSF
                                hyper_res[prm] = np.zeros((len(lutdw[lut]['meta']['wave']),len(ai))) + np.nan
                                for iii in range(len(ai)):
                                    if len(xi[0]) == 1:
                                        hyper_res[prm][:,iii] = lutdw[lut]['rgi']((xi[0], lutdw[lut]['ipd'][prm],
                                                                lutdw[lut]['meta']['wave'], xi[1], xi[2], xi[3], xi[4], ai[iii])).flatten()
                                    else:
                                        hyper_res[prm][:,iii] = lutdw[lut]['rgi']((xi[0].flatten()[iii], lutdw[lut]['ipd'][prm],
                                                                lutdw[lut]['meta']['wave'], xi[1].flatten()[iii], xi[2].flatten()[iii], xi[3].flatten()[iii], xi[4].flatten()[iii], ai[iii])).flatten()

                    ## resample to current band
                    ### path reflectance
                    romix[ls] = ac.shared.rsr_convolute_nd(hyper_res[par], lutdw[lut]['meta']['wave'], rsrd['rsr'][b]['response'], rsrd['rsr'][b]['wave'], axis=0)
                    ## transmittance and spherical albedo
                    astot[ls] = ac.shared.rsr_convolute_nd(hyper_res['astot'], lutdw[lut]['meta']['wave'], rsrd['rsr'][b]['response'], rsrd['rsr'][b]['wave'], axis=0)
                    dutott[ls] = ac.shared.rsr_convolute_nd(hyper_res['dutott'], lutdw[lut]['meta']['wave'], rsrd['rsr'][b]['response'], rsrd['rsr'][b]['wave'], axis=0)
                    ## total transmittance
                    if (setu['dsf_residual_glint_correction']) & (setu['dsf_residual_glint_correction_method']=='default'):
                        ttot_all[b][ls] = ac.shared.rsr_convolute_nd(hyper_res['ttot'], lutdw[lut]['meta']['wave'], rsrd['rsr'][b]['response'], rsrd['rsr'][b]['wave'], axis=0)
                    if (setu['output_ed']): dtott[ls] = ac.shared.rsr_convolute_nd(hyper_res['dtott'], lutdw[lut]['meta']['wave'], rsrd['rsr'][b]['response'], rsrd['rsr'][b]['wave'], axis=0)
                else:
                    ## path reflectance
                    romix[ls] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd'][par], xi[1], xi[2], xi[3], xi[4], ai))
                    ## transmittance and spherical albedo
                    astot[ls] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd']['astot'], xi[1], xi[2], xi[3], xi[4], ai))
                    dutott[ls] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd']['dutott'], xi[1], xi[2], xi[3], xi[4], ai))
                    ## total transmittance
                    if (setu['dsf_residual_glint_correction']) & (setu['dsf_residual_glint_correction_method']=='default'):
                        ttot_all[b][ls] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd']['ttot'], xi[1], xi[2], xi[3], xi[4], ai))
                    if (setu['output_ed']): dtott[ls] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd']['dtott'], xi[1], xi[2], xi[3], xi[4], ai))
                del ls, ai, xi

            ## interpolate tiled processing to full scene
            if setu['dsf_aot_estimate'] == 'tiled':
                if setu['verbosity'] > 1: print('Interpolating tiles')
                romix = ac.shared.tiles_interp(romix, xnew, ynew, target_mask=(valid_mask if setu['slicing'] else None), \
                target_mask_full=True, smooth=setu['dsf_tile_smoothing'], kern_size=setu['dsf_tile_smoothing_kernel_size'], method=setu['dsf_tile_interp_method'])
                astot = ac.shared.tiles_interp(astot, xnew, ynew, target_mask=(valid_mask if setu['slicing'] else None), \
                target_mask_full=True, smooth=setu['dsf_tile_smoothing'], kern_size=setu['dsf_tile_smoothing_kernel_size'], method=setu['dsf_tile_interp_method'])
                dutott = ac.shared.tiles_interp(dutott, xnew, ynew, target_mask=(valid_mask if setu['slicing'] else None), \
                target_mask_full=True, smooth=setu['dsf_tile_smoothing'], kern_size=setu['dsf_tile_smoothing_kernel_size'], method=setu['dsf_tile_interp_method'])
                if (setu['output_ed']):
                    dtott = ac.shared.tiles_interp(dtott, xnew, ynew, target_mask=(valid_mask if setu['slicing'] else None), \
                                                    target_mask_full=True, smooth=setu['dsf_tile_smoothing'], kern_size=setu['dsf_tile_smoothing_kernel_size'], method=setu['dsf_tile_interp_method'])

            ## create full scene parameters for segmented processing
            if setu['dsf_aot_estimate'] == 'segmented':
                romix_ = romix * 1.0
                astot_ = astot * 1.0
                dutott_ = dutott * 1.0
                romix = np.zeros(gem.gatts['data_dimensions']) + np.nan
                astot = np.zeros(gem.gatts['data_dimensions']) + np.nan
                dutott = np.zeros(gem.gatts['data_dimensions']) + np.nan
                if (setu['output_ed']):
                    dtott_ = dtott * 1.0
                    dtott = np.zeros(gem.gatts['data_dimensions']) + np.nan
                for sidx, segment in enumerate(segment_data):
                    romix[segment_data[segment]['sub']] = romix_[sidx]
                    astot[segment_data[segment]['sub']] = astot_[sidx]
                    dutott[segment_data[segment]['sub']] = dutott_[sidx]
                    if (setu['output_ed']): dtott[segment_data[segment]['sub']] = dtott_[sidx]
                del romix_, astot_, dutott_
                if (setu['output_ed']):  del dtott_

            ## write ac parameters
            if setu['dsf_write_tiled_parameters']:
                if len(np.atleast_1d(romix)>1):
                    if romix.shape == cur_data.shape:
                        gemo.write('romix_{}'.format(gem.bands[b]['wave_name']), romix)
                    else:
                        ds_att['romix'] = romix[0]
                if len(np.atleast_1d(astot)>1):
                    if astot.shape == cur_data.shape:
                        gemo.write('astot_{}'.format(gem.bands[b]['wave_name']), astot)
                    else:
                        ds_att['astot'] = astot[0]
                if len(np.atleast_1d(dutott)>1):
                    if dutott.shape == cur_data.shape:
                        gemo.write('dutott_{}'.format(gem.bands[b]['wave_name']), dutott)
                    else:
                        ds_att['dutott'] = dutott[0]

            ## do atmospheric correction
            rhot_noatm = (cur_data / gem.bands[b]['tt_gas']) - romix
            del romix
            cur_data = (rhot_noatm) / (dutott + astot*rhot_noatm)

            ## perform FFSS correction after a/c
            if (setu['dsf_interface_reflectance']) & (setu['dsf_interface_option']  == 'ffss_boa'):
                xi = (gem.data_mem['sza'+gk], gem.data_mem['vza'+gk_vza], np.abs(180-gem.data_mem['raa'+gk_raa]), aot_lut, aot_sel)
                vza = np.atleast_1d(gem.data_mem['vza'+gk_vza])
                vza[vza == 0.0] = 0.001
                rhof = ac.ac.sky_refl(np.radians(vza), n_w=1.34)
                del vza
                rsky = rhof * rgi_rsky[b](xi)
                del xi, rhof
                if (setu['dsf_aot_estimate'] == 'tiled'):
                    if setu['verbosity'] > 1: print('Interpolating tiles for rsky')
                    rsky = ac.shared.tiles_interp(rsky, xnew, ynew, target_mask=(valid_mask if setu['slicing'] else None), \
                                target_mask_full=True, smooth=setu['dsf_tile_smoothing'], kern_size=setu['dsf_tile_smoothing_kernel_size'], method=setu['dsf_tile_interp_method'])

                cur_data -= rsky
                del rsky
            ## end FFSS

            ## compute at surface Ed
            if setu['output_ed']:
                if 'sza' in gem.datasets:
                    sza = gem.data('sza')
                else:
                    sza = gem.gatts['sza']
                se_distance = ac.shared.sun_position(gem.gatts['isodate'], 0, 0)['distance']
                Ed = 1 / (1 - cur_data * astot)
                Ed *= (gem.bands[b]['F0']) * se_distance**2 * np.cos(np.radians(sza)) * gem.bands[b]['td_gas'] * dtott
                del sza
            del astot, dutott, rhot_noatm

        ## exponential
        elif (ac_opt == 'exp'):
            ## get Rayleigh correction
            rorayl_cur = lutdw[exp_lut]['rgi'][b]((xi[0], lutdw[exp_lut]['ipd'][par], xi[1], xi[2], xi[3], xi[4], 0.001))
            dutotr_cur = lutdw[exp_lut]['rgi'][b]((xi[0], lutdw[exp_lut]['ipd']['dutott'], xi[1], xi[2], xi[3], xi[4], 0.001))

            ## get epsilon in current band
            delta = (long_wv-gem.bands[b]['wave_nm'])/(long_wv-short_wv)
            eps_cur = np.power(epsilon, delta)
            rhoam_cur = rhoam * eps_cur

            ## add results to band
            if exp_fixed_epsilon: ds_att['epsilon'] = eps_cur
            if exp_fixed_rhoam: ds_att['rhoam'] = rhoam_cur

            cur_data = (cur_data - rorayl_cur - rhoam_cur) / (dutotr_cur)
            cur_data[mask] = np.nan
        ## end exponential

        ## write rhorc
        if (setu['output_rhorc']):
            ## read TOA
            cur_rhorc, cur_att = gem.data(dsi, attributes=True)

            ## compute Rayleigh parameters for DSF
            if (ac_opt == 'dsf'):
                ## no subset
                xi = [gem.data_mem['pressure'+gk],
                      gem.data_mem['raa'+gk_raa],
                      gem.data_mem['vza'+gk_vza],
                      gem.data_mem['sza'+gk],
                      gem.data_mem['wind'+gk]]

                ## get Rayleigh parameters
                if hyper:
                    rorayl_hyper = lutdw[luts[0]]['rgi']((xi[0], lutdw[luts[0]]['ipd'][par],
                                        lutdw[luts[0]]['meta']['wave'], xi[1], xi[2], xi[3], xi[4], 0.001)).flatten()
                    dutotr_hyper = lutdw[luts[0]]['rgi']((xi[0], lutdw[luts[0]]['ipd']['dutott'],
                                        lutdw[luts[0]]['meta']['wave'], xi[1], xi[2], xi[3], xi[4], 0.001)).flatten()
                    rorayl_cur = ac.shared.rsr_convolute_nd(rorayl_hyper, lutdw[luts[0]]['meta']['wave'], rsrd['rsr'][b]['response'], rsrd['rsr'][b]['wave'], axis=0)
                    dutotr_cur = ac.shared.rsr_convolute_nd(dutotr_hyper, lutdw[luts[0]]['meta']['wave'], rsrd['rsr'][b]['response'], rsrd['rsr'][b]['wave'], axis=0)
                    del rorayl_hyper, dutotr_hyper
                else:
                    rorayl_cur = lutdw[luts[0]]['rgi'][b]((xi[0], lutdw[luts[0]]['ipd'][par], xi[1], xi[2], xi[3], xi[4], 0.001))
                    dutotr_cur = lutdw[luts[0]]['rgi'][b]((xi[0], lutdw[luts[0]]['ipd']['dutott'], xi[1], xi[2], xi[3], xi[4], 0.001))
                del xi

            ## create full scene parameters for segmented processing
            if setu['dsf_aot_estimate'] == 'segmented':
                rorayl_ = rorayl_cur * 1.0
                dutotr_ = dutotr_cur * 1.0
                rorayl_cur = np.zeros(gem.gatts['data_dimensions']) + np.nan
                dutotr_cur = np.zeros(gem.gatts['data_dimensions']) + np.nan
                for sidx, segment in enumerate(segment_data):
                    rorayl_cur[segment_data[segment]['sub']] = rorayl_[sidx]
                    dutotr_cur[segment_data[segment]['sub']] = dutotr_[sidx]
                del rorayl_, dutotr_

            ## create full scene parameters for tiled processing
            if (setu['dsf_aot_estimate'] == 'tiled') & (use_revlut):
                if setu['verbosity'] > 1: print('Interpolating tiles for rhorc')
                rorayl_cur = ac.shared.tiles_interp(rorayl_cur, xnew, ynew, target_mask=(valid_mask if setu['slicing'] else None), \
                            target_mask_full=True, smooth=setu['dsf_tile_smoothing'], kern_size=setu['dsf_tile_smoothing_kernel_size'], method=setu['dsf_tile_interp_method'])
                dutotr_cur = ac.shared.tiles_interp(dutotr_cur, xnew, ynew, target_mask=(valid_mask if setu['slicing'] else None), \
                            target_mask_full=True, smooth=setu['dsf_tile_smoothing'], kern_size=setu['dsf_tile_smoothing_kernel_size'], method=setu['dsf_tile_interp_method'])

            ## write ac parameters
            if setu['dsf_write_tiled_parameters']:
                if len(np.atleast_1d(rorayl_cur)>1):
                    if rorayl_cur.shape == cur_data.shape:
                        gemo.write('rorayl_{}'.format(gem.bands[b]['wave_name']), rorayl_cur)
                    else:
                        ds_att['rorayl'] = rorayl_cur[0]
                if len(np.atleast_1d(dutotr_cur)>1):
                    if dutotr_cur.shape == cur_data.shape:
                        gemo.write('dutotr_{}'.format(gem.bands[b]['wave_name']), dutotr_cur)
                    else:
                        ds_att['dutotr'] = dutotr_cur[0]

            cur_rhorc = (cur_rhorc - rorayl_cur) / (dutotr_cur)
            gemo.write(dso.replace('rhos_', 'rhorc_'), cur_rhorc, ds_att = ds_att)
            del cur_rhorc, rorayl_cur, dutotr_cur

        if ac_opt == 'dsf' and setu['slicing']:
            del valid_mask

        ## write rhos
        gemo.write(dso, cur_data, ds_att = ds_att)
        del cur_data

        ## write Ed data
        if setu['output_ed']:
            gemo.write(dso.replace('rhos_', 'Ed_'), Ed, ds_att = ds_att)
            del Ed

        if setu['verbosity'] > 1: print('{}/B{} took {:.1f}s ({})'.format(sensor_lut, b, time.time()-t0, 'RevLUT' if use_revlut else 'StdLUT'))

    ## update outputfile dataset info
    gemo.datasets_read()

    ## glint correction
    if (ac_opt == 'dsf') & (setu['dsf_residual_glint_correction']) & (setu['dsf_residual_glint_correction_method']=='default'):
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
            sza = gem.data('sza') * dtor
            vza = gem.data('vza') * dtor
            raa = gem.data('raa') * dtor

            ## flatten 1 element arrays
            if sza.shape == (1,1): sza = sza.flatten()
            if vza.shape == (1,1): vza = vza.flatten()
            if raa.shape == (1,1): raa = raa.flatten()

            muv = np.cos(vza)
            mus = np.cos(sza)
            cos2omega = mus*muv + np.sin(sza)*np.sin(vza)*np.cos(raa)
            del sza, vza, raa

            omega = np.arccos(cos2omega)/2
            del cos2omega

            ## read and resample refractive index
            refri = ac.ac.refri()
            refri_sen = ac.shared.rsr_convolute_dict(refri['wave']/1000, refri['n'], rsrd['rsr'])

            ## compute fresnel reflectance for the reference bands
            Rf_sen = {}
            for b in [gc_swir1_b, gc_swir2_b, gc_user_b]:
                if b is None: continue
                Rf_sen[b] = ac.ac.sky_refl(omega, n_w=refri_sen[b])

            ## compute where to apply the glint correction
            ## sub_gc has the idx for non masked data with rhos_ref below the masking threshold
            gc_mask_data = gemo.data(gc_mask)

            if gc_mask_data is None: ## reference rhos dataset can be missing for night time images (tgas computation goes negative)
                print('No glint mask could be determined.')
            else:
                sub_gc = np.where(np.isfinite(gc_mask_data) & \
                                  (gc_mask_data<=setu['glint_mask_rhos_threshold']))
                del gc_mask_data

                ## get reference bands transmittance
                for ib, b in enumerate(gemo.bands):
                    rhos_ds = gemo.bands[b]['rhos_ds']
                    if rhos_ds not in [gc_swir1, gc_swir2, gc_user]: continue
                    if rhos_ds not in gemo.datasets: continue

                    ## two way direct transmittance
                    if setu['dsf_aot_estimate'] == 'tiled':
                        if setu['slicing']:
                            ## load rhos dataset
                            cur_data = gemo.data(rhos_ds)
                            valid_mask = np.isfinite(cur_data)
                            del cur_data
                        ttot_all_b = ac.shared.tiles_interp(ttot_all[b], xnew, ynew, target_mask=(valid_mask if setu['slicing'] else None), \
                        target_mask_full=True, smooth=setu['dsf_tile_smoothing'], kern_size=setu['dsf_tile_smoothing_kernel_size'], method=setu['dsf_tile_interp_method'])
                    elif setu['dsf_aot_estimate'] == 'segmented':
                        ttot_all_ = ttot_all[b] * 1.0
                        ttot_all_b = np.zeros(gem.gatts['data_dimensions']) + np.nan
                        for sidx, segment in enumerate(segment_data):
                            ttot_all_b[segment_data[segment]['sub']] = ttot_all_[sidx]
                        del ttot_all_
                    else:
                        ttot_all_b = ttot_all[b] * 1.0

                    T_cur  = np.exp(-1.*(ttot_all_b/muv)) * np.exp(-1.*(ttot_all_b/mus))
                    del ttot_all_b

                    ## subset if 2d
                    T_cur_sub = T_cur[sub_gc] if len(np.atleast_2d(T_cur)) > 1 else T_cur[0] * 1.0

                    if rhos_ds == gc_user:
                        T_USER = T_cur_sub * 1.0
                    else:
                        if rhos_ds == gc_swir1: T_SWIR1 = T_cur_sub * 1.0
                        if rhos_ds == gc_swir2: T_SWIR2 = T_cur_sub * 1.0
                    del T_cur, T_cur_sub

                ## swir band choice is made for first band
                gc_choice = False
                ## glint correction per band
                for ib, b in enumerate(gemo.bands):
                    rhos_ds = gemo.bands[b]['rhos_ds']
                    if rhos_ds not in gemo.datasets: continue
                    if b not in ttot_all: continue
                    print('Performing glint correction for band {} ({} nm)'.format(b, gemo.bands[b]['wave_name']))
                    ## load rhos dataset
                    cur_data = gemo.data(rhos_ds)

                    ## two way direct transmittance
                    if setu['dsf_aot_estimate'] == 'tiled':
                        if setu['slicing']: valid_mask = np.isfinite(cur_data)
                        ttot_all_b = ac.shared.tiles_interp(ttot_all[b], xnew, ynew, target_mask=(valid_mask if setu['slicing'] else None), \
                        target_mask_full=True, smooth=setu['dsf_tile_smoothing'], kern_size=setu['dsf_tile_smoothing_kernel_size'], method=setu['dsf_tile_interp_method'])
                    elif setu['dsf_aot_estimate'] == 'segmented':
                        ttot_all_ = ttot_all[b] * 1.0
                        ttot_all_b = np.zeros(gem.gatts['data_dimensions']) + np.nan
                        for sidx, segment in enumerate(segment_data):
                            ttot_all_b[segment_data[segment]['sub']] = ttot_all_[sidx]
                        del ttot_all_
                    else:
                        ttot_all_b = ttot_all[b] * 1.0
                    if setu['dsf_write_tiled_parameters']:
                         if len(np.atleast_1d(ttot_all_b)>1):
                             if ttot_all_b.shape == muv.shape:
                                 gemo.write('ttot_{}'.format(gem.bands[b]['wave_name']), ttot_all_b)
                             else:
                                 gemo.bands[b]['ttot_all'] = ttot_all_b[0]
                    ## end compute ttot_all_band

                    T_cur  = np.exp(-1.*(ttot_all_b/muv)) * np.exp(-1.*(ttot_all_b/mus))
                    del ttot_all_b

                    ## subset if 2d
                    T_cur_sub = T_cur[sub_gc] if len(np.atleast_2d(T_cur)) > 1 else T_cur[0] * 1.0
                    del T_cur

                    ## get current band Fresnel reflectance
                    Rf_sen_cur = ac.ac.sky_refl(omega, n_w=refri_sen[b])

                    ## get gc factors for this band
                    if gc_user is None:
                        if len(np.atleast_2d(Rf_sen[gc_swir1_b]))>1: ## if resolved angles
                            gc_SWIR1 = (T_cur_sub/T_SWIR1) * (Rf_sen_cur[sub_gc]/Rf_sen[gc_swir1_b][sub_gc])
                            gc_SWIR2 = (T_cur_sub/T_SWIR2) * (Rf_sen_cur[sub_gc]/Rf_sen[gc_swir2_b][sub_gc])
                        else:
                            gc_SWIR1 = (T_cur_sub/T_SWIR1) * (Rf_sen_cur/Rf_sen[gc_swir1_b])
                            gc_SWIR2 = (T_cur_sub/T_SWIR2) * (Rf_sen_cur/Rf_sen[gc_swir2_b])
                    else:
                        if len(np.atleast_2d(Rf_sen[gc_user_b]))>1: ## if resolved angles
                            gc_USER = (T_cur_sub/T_USER) * (Rf_sen_cur[sub_gc]/Rf_sen[gc_user_b][sub_gc])
                        else:
                            gc_USER = (T_cur_sub/T_USER) * (Rf_sen_cur/Rf_sen[gc_user_b])
                    del Rf_sen_cur, T_cur_sub

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
                            del g1_blue, g2_blue
                            rhog_ref = swir2_rhos
                            rhog_ref[use_swir1] = swir1_rhos[use_swir1]
                            del swir1_rhos, swir2_rhos
                        else:
                            rhog_ref = gemo.data(gc_user)[sub_gc]
                            ## set negatives to 0
                            rhog_ref[rhog_ref<0] = 0
                        ## write reference glint
                        if setu['glint_write_rhog_ref']:
                            tmp = np.zeros(gemo.gatts['data_dimensions'], dtype=np.float32) + np.nan
                            tmp[sub_gc] = rhog_ref
                            gemo.write('rhog_ref', tmp)
                            del tmp
                    ## end select glint correction band

                    ## calculate glint in this band
                    if gc_user is None:
                        cur_rhog = gc_SWIR2 * rhog_ref
                        try:
                            cur_rhog[use_swir1] = gc_SWIR1[use_swir1] * rhog_ref[use_swir1]
                        except:
                            cur_rhog[use_swir1] = gc_SWIR1 * rhog_ref[use_swir1]
                        del gc_SWIR1, gc_SWIR2
                    else:
                        cur_rhog = gc_USER * rhog_ref
                        del gc_USER

                    ## remove glint from rhos
                    cur_data[sub_gc]-=cur_rhog
                    gemo.write(rhos_ds, cur_data, ds_att = gem.bands[b])
                    del cur_data

                    ## write band glint
                    if setu['glint_write_rhog_all']:
                        tmp = np.zeros(gemo.gatts['data_dimensions'], dtype=np.float32) + np.nan
                        tmp[sub_gc] = cur_rhog
                        gemo.write('rhog_{}'.format(gemo.bands[b]['wave_name']), tmp, ds_att={'wavelength':gemo.bands[b]['wavelength']})
                        del tmp
                    del cur_rhog
                del sub_gc, rhog_ref
                if gc_user is not None:
                    del T_USER
                else:
                    del T_SWIR1, T_SWIR2, use_swir1
            del Rf_sen, omega, muv, mus
        if (setu['dsf_aot_estimate'] == 'tiled') & (setu['slicing']):
            del valid_mask
    ## end glint correction

    ## alternative glint correction
    if (ac_opt == 'dsf') & (setu['dsf_residual_glint_correction']) & (setu['dsf_aot_estimate'] in ['fixed', 'segmented']) &\
       (setu['dsf_residual_glint_correction_method']=='alternative'):

        ## reference aot and wind speed
        if setu['dsf_aot_estimate'] == 'fixed':
            gc_aot = max(0.1, gemo.gatts['ac_aot_550'])
            gc_wind = 20
            gc_lut = gemo.gatts['ac_model']

            raa = gem.gatts['raa'] if 'raa' in gem.gatts else np.nanmean(gem.data('raa'))
            sza = gem.gatts['sza'] if 'sza' in gem.gatts else np.nanmean(gem.data('sza'))
            vza = gem.gatts['vza'] if 'vza' in gem.gatts else np.nanmean(gem.data('vza'))

            ## flatten 1 element arrays
            if raa.shape == (1,1): raa = raa.flatten()
            if sza.shape == (1,1): sza = sza.flatten()
            if vza.shape == (1,1): vza = vza.flatten()

            ## get surface reflectance for fixed geometry
            if len(np.atleast_1d(raa)) == 1:
                if hyper:
                    surf = lutdw[gc_lut]['rgi']((gem.gatts['pressure'],lutdw[gc_lut]['ipd']['rsky_s'],
                                                 lutdw[gc_lut]['meta']['wave'], raa, vza, sza, gc_wind, gc_aot))
                    surf_res = ac.shared.rsr_convolute_dict(lutdw[gc_lut]['meta']['wave'], surf, rsrd['rsr'])
                else:
                    surf_res = {b : lutdw[gc_lut]['rgi'][b]((gem.gatts['pressure'],lutdw[gc_lut]['ipd']['rsky_s'],
                                                             raa, vza, sza, gc_wind, gc_aot)) for b in lutdw[gc_lut]['rgi']}

        if setu['dsf_aot_estimate'] == 'segmented':
            for sidx, segment in enumerate(segment_data):
                gc_aot = max(0.1, aot_sel[sidx])
                gc_wind = 20
                gc_lut = luts[aot_lut[sidx][0]]

                if sidx == 0: surf_res = {}
                ## get surface reflectance for segmented geometry
                #if len(np.atleast_1d(raa)) == 1:
                if hyper:
                    surf = lutdw[gc_lut]['rgi']((gem.data_mem['pressure'+gk][sidx],lutdw[gc_lut]['ipd']['rsky_s'],
                                                 lutdw[gc_lut]['meta']['wave'],
                                                             gem.data_mem['raa'+gk_raa][sidx],gem.data_mem['vza'+gk_vza][sidx],
                                                             gem.data_mem['sza'+gk][sidx], gc_wind, gc_aot))
                    surf_res[segment] = ac.shared.rsr_convolute_dict(lutdw[gc_lut]['meta']['wave'], surf, rsrd['rsr'])
                else:
                    surf_res[segment] = {b : lutdw[gc_lut]['rgi'][b]((gem.data_mem['pressure'+gk][sidx],lutdw[gc_lut]['ipd']['rsky_s'],
                                                            gem.data_mem['raa'+gk_raa][sidx],gem.data_mem['vza'+gk_vza][sidx],
                                                            gem.data_mem['sza'+gk][sidx], gc_wind, gc_aot)) for b in lutdw[gc_lut]['rgi']}

        ## get reference surface reflectance
        gc_ref = None
        for ib, b in enumerate(gemo.bands):
            rhos_ds = gemo.bands[b]['rhos_ds']
            if rhos_ds not in gemo.datasets: continue
            if (gemo.bands[b]['wavelength'] < setu['dsf_residual_glint_wave_range'][0]) |\
               (gemo.bands[b]['wavelength'] > setu['dsf_residual_glint_wave_range'][1]): continue
            print('Reading reference for glint correction from band {} ({} nm)'.format(b, gemo.bands[b]['wave_name']))

            if setu['dsf_aot_estimate'] == 'fixed':
                gc_sur_cur = surf_res[b]
            if setu['dsf_aot_estimate'] == 'segmented':
                gc_sur_cur = gemo.data(rhos_ds) * np.nan
                for segment in segment_data:
                    gc_sur_cur[segment_data[segment]['sub']] = surf_res[segment][b]

            if gc_ref is None:
                gc_ref = gemo.data(rhos_ds)
                gc_sur = gc_sur_cur
            else:
                gc_ref = np.dstack((gc_ref, gemo.data(rhos_ds)))
                gc_sur = np.dstack((gc_sur, gc_sur_cur))

        if gc_ref is None:
            print('No bands found between {} and {} nm for glint correction'.format(setu['dsf_residual_glint_wave_range'][0],
                                                                                    setu['dsf_residual_glint_wave_range'][1]))
        else:
            ## compute average reference glint
            if len(gc_ref.shape) == 3:
                gc_ref[gc_ref<0] = np.nan
                gc_ref_mean = np.nanmean(gc_ref, axis=2)
                gc_ref_std = np.nanstd(gc_ref, axis=2)
                gc_ref = None

                gemo.write('glint_mean', gc_ref_mean)
                gemo.write('glint_std', gc_ref_std)
            else: ## or use single band
                gc_ref[gc_ref<0] = 0.0
                gc_ref_mean = gc_ref*1.0

            ## compute average modeled surface glint
            axis = None
            if setu['dsf_aot_estimate'] == 'segmented': axis = 2
            gc_sur_mean = np.nanmean(gc_sur, axis = axis)
            gc_sur_std = np.nanstd(gc_sur, axis = axis)

            ## get subset where to apply glint correction
            gc_sub = np.where(gc_ref_mean<setu['glint_mask_rhos_threshold'])

            ## glint correction per band
            for ib, b in enumerate(gemo.bands):
                rhos_ds = gemo.bands[b]['rhos_ds']
                if rhos_ds not in gemo.datasets: continue
                print('Performing glint correction for band {} ({} nm)'.format(b, gemo.bands[b]['wave_name']))

                ## estimate current band glint from reference glint image and ratio of interface reflectance
                if setu['dsf_aot_estimate'] == 'fixed':
                    sur = surf_res[b] * 1.0

                if setu['dsf_aot_estimate'] == 'segmented':
                    sur = gc_ref_mean * np.nan
                    for segment in segment_data:
                        sur[segment_data[segment]['sub']] = surf_res[segment][b]
                    sur = sur[gc_sub]

                if len(np.atleast_2d(gc_sur_mean)) == 1:
                    cur_rhog = gc_ref_mean[gc_sub] * (sur/gc_sur_mean)
                else:
                    cur_rhog = gc_ref_mean[gc_sub] * (sur/gc_sur_mean[gc_sub])

                ## remove glint from rhos
                cur_data = gemo.data(rhos_ds)
                cur_data[gc_sub]-=cur_rhog
                gemo.write(rhos_ds, cur_data, ds_att = gem.bands[b])

                ## write band glint
                if setu['glint_write_rhog_all']:
                    tmp = np.zeros(gemo.gatts['data_dimensions'], dtype=np.float32) + np.nan
                    tmp[gc_sub] = cur_rhog
                    gemo.write('rhog_{}'.format(gemo.bands[b]['wave_name']), tmp, ds_att={'wavelength':gemo.bands[b]['wavelength']})
                    tmp = None
                cur_rhog = None
    ## end alternative glint correction

    ## compute contrabands
    if setu['compute_contrabands']: ac.parameters.castagna.contraband(gemo)

    ## clear aot results
    aot_lut, aot_sel = None, None

    ## close inputfile
    gem.close()
    gem = None

    ## update attributes with latest version
    if output_file:
        gemo.gatts_update()
        gemo.close()

    if setu['verbosity']>0: print('Wrote {}'.format(ofile))

    if return_gem:
        return(gemo, setu)
    else:

        return(ofile, setu)
