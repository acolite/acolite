## def dsf_optimise_aot
## function to determine optimal model and aot for homogeneous processing given L1R file
## returns selected lut and aot e.g. for acolite_l2r "DSF" processing
## written by Quinten Vanhellemont, RBINS
## 2025-02-10
## modifications: 2025-03-17 (QV) moved target reading to file, compute average geometry, check if target lon and lat are given,
##                                fix for interface reflectance, added sensor/detector name, added support for hyper sensors

def optimise_aot_homogeneous(gem, quiet = True, settings = None, lutdw = None):
    import numpy as np
    import scipy.optimize
    import acolite as ac

    ### correct_band_homogeneous
    ## nested function to perform correction
    ## QV 2025-02-10 adapted from RAdCor correct_band, removed writing option
    def correct_band_homogeneous(b, aot):
        ## Load rhot data
        if not quiet: print('Loading TOA data for band {} ({})'.format(b, bands[b]['rhot_ds']))
        rho_toa, att = gem.data(bands[b]['rhot_ds'], attributes = True)
        rho_toa_mask = np.isnan(rho_toa)

        ## Get atmospheric parameters
        if not hyper:
            if add_rsky:
                rho_a     = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd'][romix_par], raa, vza, sza, wind, aot))
                rho_a_sph = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['astot'], raa, vza, sza, wind, aot))
                T_du_tot   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['dutott'], raa, vza, sza, wind, aot))
            else:
                rho_a     = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['romix'], raa, vza, sza, aot))
                rho_a_sph = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['astot'], raa, vza, sza, aot))
                T_du_tot   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['dutott'], raa, vza, sza, aot))
        else:
            if add_rsky:
                rho_a     = lutdw[lut]['rgi']((pressure, lutdw[lut]['ipd'][romix_par], lutdw[lut]['meta']['wave'], raa, vza, sza, wind, aot))
                rho_a_sph = lutdw[lut]['rgi']((pressure, lutdw[lut]['ipd']['astot'], lutdw[lut]['meta']['wave'],raa, vza, sza, wind, aot))
                T_du_tot   = lutdw[lut]['rgi']((pressure, lutdw[lut]['ipd']['dutott'], lutdw[lut]['meta']['wave'],raa, vza, sza, wind, aot))
            else:
                rho_a     = lutdw[lut]['rgi']((pressure, lutdw[lut]['ipd']['romix'], lutdw[lut]['meta']['wave'],raa, vza, sza, aot))
                rho_a_sph = lutdw[lut]['rgi']((pressure, lutdw[lut]['ipd']['astot'], lutdw[lut]['meta']['wave'],raa, vza, sza, aot))
                T_du_tot   = lutdw[lut]['rgi']((pressure, lutdw[lut]['ipd']['dutott'], lutdw[lut]['meta']['wave'],raa, vza, sza, aot))

            ## resample
            rho_a = ac.shared.rsr_convolute_nd(rho_a, lutdw[lut]['meta']['wave'], rsrd[sensor_lut]['rsr'][b]['response'], rsrd[sensor_lut]['rsr'][b]['wave'], axis = 0)
            rho_a_sph = ac.shared.rsr_convolute_nd(rho_a_sph, lutdw[lut]['meta']['wave'], rsrd[sensor_lut]['rsr'][b]['response'], rsrd[sensor_lut]['rsr'][b]['wave'], axis = 0)
            T_du_tot = ac.shared.rsr_convolute_nd(T_du_tot, lutdw[lut]['meta']['wave'], rsrd[sensor_lut]['rsr'][b]['response'], rsrd[sensor_lut]['rsr'][b]['wave'], axis = 0)

        ## Gas transmittance correction
        rho_toa /= bands[b]['tt_gas']

        ## homogeneous surface correction
        rho_s_homo = rho_toa - rho_a
        rho_s_homo /= (T_du_tot + rho_s_homo * rho_a_sph)
        return(rho_s_homo)
    ### end of correct_band_homogeneous function

    ## aot optimisation function
    ## QV 2025-02-10 adapted from RAdCor version
    def opt_aot(aot, return_res = False, correct_band = correct_band_homogeneous):
        print('    {} {}'.format(lut, aot))
        ## make temporary gem
        gemt = ac.gem.gem(None, new = False)
        gemt.gatts = {k: gem.gatts[k] for k in gem.gatts}
        gemt.gatts['uoz'] = uoz
        gemt.gatts['uwv'] = uwv
        gemt.gatts['pressure'] = pressure
        gemt.gatts['wind'] = wind

        ## run through bands running correct band for current model and aot
        opt_bands_ = [b for b in opt_bands]
        if setu['dsf_residual_glint_correction']: opt_bands_ = [b for b in bands_] ##
        for bi, b in enumerate(opt_bands_):
            ds = bands[b]['rhos_ds']
            rhos_ = correct_band(b, aot)
            gemt.data_mem[ds] = rhos_

        ## perform glint correction
        if setu['dsf_residual_glint_correction']:
            ## set current aot and model
            gemt.gatts['ac_aot_550'] = aot
            gemt.gatts['ac_model'] = am
            ## add data_mem datasets
            gemt.datasets = [k for k in gemt.data_mem.keys()]
            ## run glint correction
            ret = ac.glint.default(gemt, settings = setu, new_file = False, write = False, lutdw = lutdw)

        ## extract result
        res_rhos = np.zeros(len(opt_bands))
        for bi, b in enumerate(opt_bands):
            ds = bands[b]['rhos_ds']
            if setu['optimise_target_type'] == 'pixel':
                res_rhos[bi] = gemt.data_mem[ds][i,j]
            else:
                res_rhos[bi] = np.nanmean(gemt.data_mem[ds][optimise_mask])
        del gemt ## delete temporary gem

        if return_res:
            return(res_rhos)
        else:
            if optimise_aot_cost == 'RMSD':
                return(np.nanmean((opt_rhos-res_rhos)**2)**0.5)
            elif optimise_aot_cost == 'MAPD':
                mn = 1.0 * opt_rhos ## MAPD
                mn2 = np.abs(opt_rhos-res_rhos) / mn
                return(np.nansum(mn2) / len(opt_bands))
        ## end aot optimisation function

    ## determine which dataset was passed
    opened = False
    if type(gem) is str:
        gem = ac.gem.gem(gem)
        opened = True
    else:
        if gem.file is not None:
            gem.setup() ## update dataset info
    gemf = gem.file

    ## get sensor
    sensor = gem.gatts['sensor']

    ## get run/user/sensor settings
    setu = ac.acolite.settings.merge(sensor = gem.gatts['sensor'], settings = settings)

    if (setu['optimise_target_lon'] is None):
        print('Please provide optimise_target_lon for optimisation.')
        if opened: gem = None
        return

    if (setu['optimise_target_lat'] is None):
        print('Please provide optimise_target_lat for optimisation.')
        if opened: gem = None
        return

    ## determine which LUT to load
    if (setu['dsf_interface_reflectance']):
        if (setu['dsf_interface_option'] == 'default'):
            romix_par = 'romix+rsky_t'
        elif (setu['dsf_interface_option']  == '6sv'):
            romix_par = 'romix+rsurf'
    else:
        romix_par = 'romix'
    add_rsky = romix_par != 'romix'

    ## get sensor lut
    if setu['rsr_version'] is not None:
        sensor_lut = '{}_{}'.format(sensor, setu['rsr_version'])
    else:
        sensor_lut = '{}'.format(sensor)

    ## Load RSR dict
    hyper = False
    if sensor in ac.config['hyper_sensors']: ## to add PACE/OCI SWIR RSR
        hyper = True
        rsr = ac.shared.rsr_hyper(gem.gatts['band_waves'], gem.gatts['band_widths'], step=0.1)
        rsrd = ac.shared.rsr_dict(rsrd = {sensor_lut : {'rsr' : rsr}})
    else:
        rsrd = ac.shared.rsr_dict(sensor = sensor_lut)

    ## get ancillary data
    anc = {}
    if setu['ancillary_data']:
        clon = np.nanmedian(gem.data('lon'))
        clat = np.nanmedian(gem.data('lat'))
        print('Getting ancillary data for {:.2f}N {:.2f}E at {}'.format(clat, clon, gem.gatts['isodate']))
        anc  = ac.ac.ancillary.get(gem.gatts['isodate'], clon, clat)

    ## gas concentration and pressure
    if 'uoz' not in gem.gatts:
        uoz = setu['uoz_default'] * 1.0
        if setu['ancillary_data'] & ('uoz' in anc): uoz = anc['uoz'] * 1.0 ## new anc format is already converted QV 2024-02-12
    else: uoz = gem.gatts['uoz']

    if 'uwv' not in gem.gatts:
        uwv = setu['uwv_default'] * 1.0
        if setu['ancillary_data'] & ('uwv' in anc): uwv = anc['uwv'] * 1.0 ## new anc format is already converted QV 2024-02-12
    else: uwv = gem.gatts['uwv']

    if 'pressure' not in gem.gatts:
        pressure = setu['pressure'] * 1.0
        if setu['ancillary_data'] & ('pressure' in anc): pressure = anc['pressure'] * 1.0 ## new anc format is already converted QV 2024-02-12
    else: pressure = gem.gatts['pressure']

    if 'wind' not in gem.gatts:
        if setu['wind'] is None:
            wind = setu['wind_default'] * 1.0
        else:
            wind = setu['wind'] * 1.0
        if setu['ancillary_data'] & ('wind' in anc): wind = anc['wind'] * 1.0
    else: wind = gem.gatts['wind']

    ## elevation derived pressure
    elevation = None
    if setu['elevation'] is not None:
        elevation = float(setu['elevation'])

    if setu['dem_pressure']:
        if setu['verbosity'] > 1: print('Extracting {} DEM data'.format(setu['dem_source']))
        if (('lon' in gem.datasets) & ('lat' in gem.datasets)):
            elevation = ac.dem.dem_lonlat(gem.data('lon'), gem.data('lat'), source = setu['dem_source'])
        elif (('lat' in gem.gatts) & ('lon' in gem.gatts)):
            elevation = ac.dem.dem_lonlat(gem.gatts('lon'), gem.gatts('lat'), source = setu['dem_source'])
        else:
            if setu['verbosity'] > 1: print('No latitude and longitude in file for extracting DEM data')

    if elevation is not None:
        median_elevation = np.nanmedian(elevation)
        pressure = ac.ac.pressure_elevation(median_elevation)
        print('Using {:.2f} hPa pressure derived from {:.1f} m elevation'.format(pressure, median_elevation))
    ## end elevation derived pressure

    print(', '.join(['{}'.format(v) for v in ['uoz', 'uwv', 'pressure', 'wind']]) + ': ' +\
          ', '.join(['{:.2f}'.format(v) for v in [uoz, uwv, pressure, wind]]))

    ## get average geometry
    sza = gem.gatts['sza'] if 'sza' in gem.gatts else np.nanmean(gem.data('sza'))
    vza = gem.gatts['vza'] if 'vza' in gem.gatts else np.nanmean(gem.data('vza'))
    raa = gem.gatts['raa'] if 'raa' in gem.gatts else np.nanmean(gem.data('raa'))

    ## available rhot datasets
    rhot_datasets = [ds for ds in gem.datasets if ds.startswith('rhot_')]

    ## extract point
    ## ext = ac.shared.nc_extract_point(ncf, setu['optimise_target_lon'],  setu['optimise_target_lat'], box_size = setu['optimise_target_size'])

    ## if we have resolved geometry we can use the pixel values
    ## eg from ext above

    ## get gas transmittance
    tg_dict = ac.ac.gas_transmittance(sza, vza, uoz = uoz, uwv = uwv, rsr = rsrd[sensor_lut]['rsr'])

    ## Bands dictionary
    bands = {}
    for bi, b in enumerate(rsrd[sensor_lut]['rsr_bands']):
        if b not in bands:
            bands[b] = {k:rsrd[sensor_lut][k][b] for k in ['wave_mu', 'wave_nm', 'wave_name'] if b in rsrd[sensor_lut][k]}
            bands[b]['rhot_ds'] = 'rhot_{}'.format(bands[b]['wave_name'])
            bands[b]['rhos_ds'] = 'rhos_{}'.format(bands[b]['wave_name'])

            if setu['add_band_name']:
                bands[b]['rhot_ds'] = 'rhot_{}_{}'.format(b, bands[b]['wave_name'])
                bands[b]['rhos_ds'] = 'rhos_{}_{}'.format(b, bands[b]['wave_name'])
            if setu['add_detector_name']:
                dsname = rhot_datasets[bi][5:]
                bands[b]['rhot_ds'] = 'rhot_{}'.format(dsname)
                bands[b]['rhos_ds'] = 'rhos_{}'.format(dsname)

            if bands[b]['rhot_ds'] not in gem.datasets:
                print('{} dataset for band {} not in inputfile.'.format(bands[b]['rhot_ds'], b))
                continue

            for k in tg_dict:
                if k not in ['wave']: bands[b][k] = tg_dict[k][b]
            bands[b]['wavelength'] = bands[b]['wave_nm']

            bands[b]['dsf_use_band'] = True
            if bands[b]['tt_gas'] < setu['min_tgas_rho']: bands[b]['dsf_use_band'] = False
            if bands[b]['tt_gas'] < setu['min_tgas_aot']: bands[b]['dsf_use_band'] = False

    ## List of bands to be used
    bands_ = [b for b in bands if bands[b]['dsf_use_band']]
    if len(bands_) == 0:
        print('No band data found in inputfile.')
        if opened: gem = None
        return
    else:
        print('Band list for sensor {}: {}'.format(sensor, ', '.join([str(b) for b in bands_])))

    ## if optimise_target_rhos_file is given, read the data, and resample to the sensor RSR
    if (setu['optimise_target_rhos_file'] is not None):
        ## read file for optimisation
        ret = ac.ac.optimise_read(setu)
        if ret is None:
            if opened: gem = None
            return
        else:
            wave_data, rhos_data = ret

        ## convolute to sensor bands
        rhos_bands = ac.shared.rsr_convolute_dict(wave_data/1000, rhos_data,  rsrd[sensor_lut]['rsr'], fill_value = np.nan)
        opt_rhos = np.asarray([rhos_bands[b] for b in bands_])
        print('Setting optimise_target_rhos to: {}'.format(', '.join(['{:.5f}'.format(v) for v in opt_rhos])))
        setu['optimise_target_rhos'] = opt_rhos

    ## check if number of target rhos corresponds to the number of considered bands
    if (setu['optimise_target_rhos'] is None):
        print('The setting optimise_target_rhos is None, provide optimise_target_rhos for each considered band: {}'.format(bands_))
        print('Set missing bands (e.g. SWIR) to NaN to be ignored, or to 0 to take them into account in the fit')
        if opened: gem = None
        return

    ## check if number of target rhos corresponds to the number of considered bands
    if len(setu['optimise_target_rhos']) != len(bands_):
        print('The number of items in optimise_target_rhos ({}) does not match the number of considered bands ({}).'.format(len(setu['optimise_target_rhos']), len(bands_)))
        print('Provide optimise_target_rhos for each considered band: {}'.format(bands_))
        print('Set missing bands (e.g. SWIR) to NaN to be ignored, or to 0 to take them into account in the fit')
        if opened: gem = None
        return

    ## don't provide all NaNs!
    if not any(np.isfinite(np.asarray(setu['optimise_target_rhos'], dtype=float))):
        print('Zero finite items in optimise_target_rhos: {}'.format(', '.join([str(v) for v in setu['optimise_target_rhos']])))
        if opened: gem = None
        return

    ## get target bands and rhos
    opt_rhos = np.asarray(setu['optimise_target_rhos'], dtype = np.float32)
    print('The number of items in optimise_target_rhos ({}) matches the number of considered bands ({}).'.format(len(setu['optimise_target_rhos']), len(bands_)))
    print('The provided optimise_target_rhos for each considered band:')
    print(', '.join(['{}: {}'.format(b, opt_rhos[bi]) for bi, b in enumerate(bands_)]))
    opt_sub = np.where(np.isfinite(opt_rhos))
    opt_rhos = opt_rhos[opt_sub]
    opt_bands = [bands_[s] for s in opt_sub[0]]
    print('Optimising aot to bands: {}'.format(', '.join([str(v) for v in opt_bands])))
    print('Optimising aot to rhos: {}'.format(', '.join([str(v) for v in opt_rhos])))

    ## get target pixel
    st_lat = setu['optimise_target_lat']
    st_lon = setu['optimise_target_lon']
    lon = gem.data('lon')
    lat = gem.data('lat')
    latrange = np.nanpercentile(lat, (0,100))
    lonrange = np.nanpercentile(lon, (0,100))
    if (st_lat < latrange[0]) | (st_lat > latrange[1]) | (st_lon < lonrange[0]) | (st_lon > lonrange[1]):
        print('Point {}N {}E not in scene {}'.format(st_lat, st_lon, ncf))
        if opened: gem = None
        return

    ## load LUTs
    if lutdw is None:
        lutdw = ac.aerlut.import_luts(add_rsky = add_rsky, par = romix_par,
                                      sensor = sensor_lut if not hyper else None,
                                      rsky_lut = setu['dsf_interface_lut'], base_luts = setu['luts'], pressures = setu['luts_pressures'],
                                      reduce_dimensions = False)
    luts = list(lutdw.keys())

    ## Find pixel
    tmp = ((lon - st_lon)**2 + (lat - st_lat)**2)**0.5
    i, j = np.where(tmp == np.nanmin(tmp))
    i = i[0]
    j = j[0]
    print('Optimising aot to lat, lon: {}, {}'.format(st_lat, st_lon))
    if setu['optimise_target_type'] == 'pixel':
        print('Optimising aot to pixel: {}, {}'.format(i, j))
        print('Pixel coordinates: {}, {}'.format(lat[i,j], lon[i,j]))
    else:
        target_dim = setu['optimise_target_size']
        if setu['optimise_target_units'][0].lower() == 'm': target_dim /= resolution

        ## compute target mask for box/radius matching
        if setu['optimise_target_type'] == 'box': ## square mask
            hbox = int(target_dim/2)
            print('Optimising aot to {} pixel {} with centre pixel: {}, {}'.format(hbox*2, setu['optimise_target_type'], i, j))
            optimise_target_mask = np.zeros(lon.shape, dtype = bool)
            optimise_target_mask[i-hbox:i+hbox+1, j-hbox:j+hbox+1] = True
        elif setu['optimise_target_type'] == 'circle': ## circular mask
            radius_pix = int(target_dim)
            print('Optimising aot to {} with {} pixel radius and centre pixel: {}, {}'.format(setu['optimise_target_type'],radius_pix, i, j))
            y = np.arange(0, lon.shape[0])
            x = np.arange(0, lon.shape[1])
            optimise_target_mask = (x[None, :] - j) ** 2 + (y[:, None] - i) ** 2 < radius_pix**2
            print('Centre pixel coordinates: {}, {}'.format(lat[i,j], lon[i,j]))

        ## add target mask with default acolite masking
        if setu['optimise_target_mask']:
            flags = ac.acolite.acolite_flags(gem, create_flags_dataset = False, write_flags_dataset = False,
                                                    return_flags_dataset = True)
            optimise_target_mask[np.where(flags != 0)] = False ## remove pixels from target mask if non-zero flags

        ## compute optimisation mask
        ## add additional masking before this step
        optimise_mask = np.where(optimise_target_mask)
        print('Subset for optimisation is {} valid target pixels.'.format(len(optimise_mask[0])))
        if len(optimise_mask[0]) == 0:
            print('Optimisation not possible with 0 valid target pixels. Please change masking settings.')
            if opened: gem = None
            return

    ## cost function for aot estimate
    optimise_aot_cost = setu['optimise_aot_cost'].upper()
    print('Optimising aot with cost function: {}'.format(optimise_aot_cost))

    ## Fit AOT(550) for each model
    aer_models = ["C", "M"]
    aer_nm = {'C': 'MOD1', 'M': 'MOD2'}

    print('    Estimating aot550 for models {}'.format(', '.join(aer_models)))
    model_band_selection = {}
    for ai, am in enumerate(aer_models):
        lut = [lut for lut in luts if aer_nm[am] in lut]
        if len(lut) == 0:
            print('    Model {} not in run luts setting'.format(am))
            continue
        else:
            lut = lut[0]
        print('    Estimating aot550 for model {}'.format(am))
        opt = scipy.optimize.minimize_scalar(opt_aot, bounds = (0.001, 5.0), method = 'bounded', options = {"xatol":setu['optimise_tolerance']})
        model_band_selection[am] = {'aot': opt.x, 'fit': opt_aot(opt.x), 'result': opt_aot(opt.x, return_res = True), 'lut': lut}
        print('    Optimised aot550 for model {}: {:.4f}'.format(am, model_band_selection[am]['aot']))
        print('    Fit: {:.4e}'.format(model_band_selection[am]['fit']))

    best_mod = None
    for ai, am in enumerate(aer_models):
        if am not in model_band_selection: continue
        print('    Model {} fit to rhos: {:.4f}'.format(am, model_band_selection[am]['fit']))
        if best_mod is None:
            best_mod = am
            best_fit = 1.0 * model_band_selection[am]['fit']
            best_idx = ai
            best_aot = model_band_selection[am]['aot']
        elif best_fit > model_band_selection[am]['fit']:
            best_mod = am
            best_fit = 1.0 * model_band_selection[am]['fit']
            best_idx = ai
            best_aot = model_band_selection[am]['aot']

    best_lut = [lut for lut in luts if aer_nm[best_mod] in lut][0]
    print('\nOptimised aerosol model: {} ({})'.format(best_mod, best_lut))
    print('Optimised aerosol optical thickness at 550 nm: {:.4f}'.format(best_aot))

    ## plot results
    if setu['optimise_plot']:
        import matplotlib.pyplot as plt

        opt_wave = np.asarray([bands[b]['wavelength'] for b in opt_bands])
        opt_sort = np.argsort(opt_wave)

        ## sel model for annotating plot
        sel_am = None
        for am in model_band_selection:
            if am not in model_band_selection: continue
            if sel_am is None:
                sel_am = am
            elif model_band_selection[am]['fit']<model_band_selection[sel_am]['fit']:
                sel_am = am

        ## plot result per model
        fig, ax = plt.subplots()
        plt.plot(opt_wave[opt_sort], opt_rhos[opt_sort], '.-', color='Black', label = 'target')
        if (setu['optimise_target_rhos_file'] is not None):
            plt.plot(wave_data, rhos_data, ':', color='Grey', label = 'target (hyperspectral)')
        for ai, am in enumerate(aer_models):
            if am not in model_band_selection: continue
            if am == 'M':
                col = '#1f77b4'
            if am == 'C':
                col = '#ff7f0e'
            if am == 'U':
                col = '#ff0000'
            plt.plot(opt_wave[opt_sort], model_band_selection[am]['result'][opt_sort],  '.:', color = col,
                    label = r'MOD{} $\tau_a$={:.3f} {}={:.2e} {}'.format(am,model_band_selection[am]['aot'],
                                                                              optimise_aot_cost, model_band_selection[am]['fit'],
                                                                              '(*)' if sel_am == am else ''))
        plt.legend()
        plt.title('{} {} DSF'.format(sensor.replace('_', '/'), gem.gatts['isodate'][0:19]))
        plt.xlabel('Wavelength (nm)')
        plt.ylabel(r'$\rho_{s}$ (1)')
        xlim = plt.xlim()
        plt.plot(xlim, [0,0], '--', color='Grey')
        plt.xlim(xlim)
        plt.ylim(setu['optimise_plot_range'][0], setu['optimise_plot_range'][1])
        plt.savefig('{}/{}_{}.png'.format(setu['output'], gem.gatts['oname'], 'rhos_optimised'), dpi = 300, bbox_inches = 'tight', facecolor = 'white')
        plt.close()

    return(best_lut, best_aot)
