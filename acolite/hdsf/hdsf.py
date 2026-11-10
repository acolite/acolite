## def hdsf
## runs lira processing based on DSF output file (only dsf_aot_estimate=fixed and ancillary_fixed supported)
## written by Quinten Vanhellemont, RBINS
## 2026-07-29
## modifications: 2026-07-29 (QV) added ancillary_fixed and dsf_aot_fixed option

def hdsf(ncf, output = None, settings = None):
    import os
    import numpy as np
    import acolite as ac

    ## get file attributes and datasets
    if (type(ncf) is str): gem = ac.gem.gem(ncf)
    gatts = gem.gatts
    datasets = gem.datasets

    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## get sensor specific defaults
    setd = ac.acolite.settings.parse(gatts['sensor'])
    ## set sensor default if user has not specified the setting
    for k in setd:
        if k not in ac.settings['user']: setu[k] = setd[k]
    ## end set sensor specific defaults

    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    if ('ac_aot_550' in gatts) & ('ac_model' in gatts):
        print('Running hybrid DSF')
        aer_depth = gatts['ac_aot_550']
        aer_model = gatts['ac_model'][-1]
    elif (gatts['acolite_file_type'] == 'L1R'):
        if (setu['dsf_aot_estimate'] == 'ancillary_fixed'):
            print('Running hybrid DSF with ancillary aerosol information')
            aer_anc = ac.ac.ancillary.aer.select(gem.gatts['isodate'], np.nanmean(gem.data('lon')), np.nanmean(gem.data('lat')))
            if aer_anc is None: return
            aer_lut, aer_aot, aer_ang_mean, aer_ssa = aer_anc
            print('Setting aer_depth={:.3f} (mean) and aer_model={} (mean angstrom={:.2f}) based on ancillary data'.format(np.nanmean(aer_aot), aer_lut, aer_ang_mean))
            aer_depth = np.nanmean(aer_aot)
            aer_model = aer_lut[-1]
        elif (setu['dsf_fixed_aot'] is not None) & (setu['dsf_fixed_lut'] is not None):
            aer_depth = setu['dsf_fixed_aot']
            aer_model = setu['dsf_fixed_lut'][-1]
        else:
            print('Not running hybrid DSF for L1R file {}'.format(ncf))
            print('Set dsf_aot_option=ancillary_fixed or user defined aerosol dsf_fixed_aot and dsf_fixed_lut')
            return
    else:
        print('Not running hybrid DSF: {} not supported'.format(ncf))
        return

    ## select S&F aerosol_haze types in libRadtran
    # 1 Rural type aerosols.
    # 4 Maritime type aerosols.
    # 5 Urban type aerosols.
    # 6 Tropospheric type aerosols.
    if aer_model == '1': ## Continental 6SV
        aer_model = 'aerosol_haze 1'
    elif aer_model == '2': ## Maritime 6SV
        aer_model = 'aerosol_haze 4'
    elif aer_model == '3': ## Maritime 6SV
        aer_model = 'aerosol_haze 5'
    else:
        aer_model = 'aerosol_haze 1'

    ## get geometry and run libradtran
    sza = None
    if 'sza' in gatts:
        sza = gatts['sza'] * 1.0
    elif 'sza' in datasets:
        sza = np.nanmean(ac.shared.nc_data(ncf, 'sza'))
    if sza is None:
        print('Could not determine sza')
        return

    vza = None
    if 'vza' in gatts:
        vza = gatts['vza'] * 1.0
    elif 'vza' in datasets:
        vza = np.nanmean(ac.shared.nc_data(ncf, 'vza'))
    if vza is None:
        print('Could not determine vza')
        return

    raa = None
    if 'raa' in gatts:
        raa = gatts['raa'] * 1.0
    elif 'raa' in datasets:
        raa = np.nanmean(ac.shared.nc_data(ncf, 'raa'))

    if raa is None:
        vaa = None
        if 'vaa' in gatts:
            vaa = gatts['vaa'] * 1.0
        elif 'vaa' in datasets:
            vaa = np.nanmean(ac.shared.nc_data(ncf, 'vaa'))
        if vaa is None:
            print('Could not determine vaa')
            return

        saa = None
        if 'saa' in gatts:
            saa = gatts['saa'] * 1.0
        elif 'saa' in datasets:
            saa = np.nanmean(ac.shared.nc_data(ncf, 'saa'))
        if saa is None:
            print('Could not determine saa')
            return

        raa = np.abs(vaa - saa)
        if raa > 180: raa -= 180

    ## pressure and wind speed
    pressure = setu['pressure_default'] * 1.0
    if 'pressure' in gatts: pressure = gatts['pressure'] * 1.0
    wind = setu['wind_default'] * 1.0
    if 'wind' in gatts: wind = gatts['wind'] * 1.0
    uoz = setu['uoz_default'] * 1.0
    if 'uoz' in gatts: uoz = gatts['uoz'] * 1.0
    uwv = setu['uwv_default'] * 1.0
    if 'uwv' in gatts: uwv = gatts['uwv'] * 1.0

    ## append config
    append_cfg = ['aerosol_default', aer_model, \
                  'aerosol_set_tau_at_wvl 550 {}'.format(aer_depth)]

    ## wavelength range
    lira_wavelength = '300 2500'
    lira_reptran = 'coarse'
    lira_wind = None
    lira_o3 = 0.3
    lira_h2o = 1.5
    lira_path = 'rho_path'
    #lira_path = 'rho_path_ocean'
    lira_sa = True

    # run standard simulation
    print('Running libRadtran simulation aerosol optical depth at 550 nm {:.3f} model {}'.format(aer_depth, aer_model))
    lr = ac.rtm.libradtran.acpar(sza, vza, raa, reptran = lira_reptran,  wavelength = lira_wavelength,
                                 wind = lira_wind, o3 = lira_o3, h2o = lira_h2o,
                                  append_cfg = append_cfg, quiet = True, return_data = False)

    ## get RSR
    sensor_lut = gatts['sensor']
    if setu['rsr_version'] is not None:
        sensor_lut = '{}_{}'.format(gatts['sensor'], setu['rsr_version'])
    hyper = True if gatts['sensor'] in ac.config['hyper_sensors'] else False

    ## create or load RSR
    if (hyper) & ('band_waves' in gatts) & ('band_widths' in gatts):
        ## make hyperspectral RSR
        rsr = ac.shared.rsr_hyper(gatts['band_waves'],
                                  gatts['band_widths'], step=0.1)
        ## update PACE OCI response with known RSR
        if gatts['sensor'] == 'PACE_OCI':
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
        rsrd = ac.shared.rsr_dict(rsrd={sensor_lut:{'rsr':rsr}})[sensor_lut]
        del rsr
    else:
        rsrd = ac.shared.rsr_dict(sensor_lut)[sensor_lut]

    ## resample results
    print('Resampling libRadtran simulation results to {}'.format(sensor_lut))
    lrr = {k: ac.shared.rsr_convolute_dict(lr['wavelength']/1000, lr[k], rsrd['rsr']) for k in lr}

    ## get gas transmittance  - only used for determining whether to compute rhos
    tg_dict = ac.ac.gas_transmittance(sza, vza, uoz = uoz, uwv = uwv, rsr = rsrd['rsr'])

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
                gem.bands[b]['rhot_ds'] = 'rhot_{}_{}'.format(gem.gatts['band_detectors'][bi], gem.bands[b]['wave_name'])
                gem.bands[b]['rhos_ds'] = 'rhos_{}_{}'.format(gem.gatts['band_detectors'][bi], gem.bands[b]['wave_name'])
            if setu['output_ed']:
                gem.bands[b]['F0'] = f0_b[b]
                gem.bands[b]['td_gas'] = tdg_b[b]
            for k in tg_dict:
                if k not in ['wave']: gem.bands[b][k] = tg_dict[k][b]
            gem.bands[b]['wavelength'] = gem.bands[b]['wave_nm']

    ## create output file
    if output is None: output = setu['output']
    gemf = gem.file
    bn = os.path.basename(gemf)
    dn = os.path.dirname(gemf)
    if '_L2R' in bn:
        oname = bn.replace('_L2R', '_L2R_hDSF')
    elif '_L1R' in bn:
        oname = bn.replace('_L1R', '_L2R_hDSF')
    if output is None:
        ofile = '{}/{}'.format(dn, oname)
    else:
        ofile = '{}/{}'.format(output, oname)

    ## set up output gem
    gemo = ac.gem.gem(ofile, new = True)
    gemo.gatts = {k: gem.gatts[k] for k in gem.gatts}
    gemo.nc_projection = gem.nc_projection
    gemo.gatts['acolite_file_type'] = 'L2R'
    gemo.gatts['ofile'] = ofile

    gemo.bands = gem.bands
    gemo.verbosity = setu['verbosity']
    gemo.gatts['acolite_version'] = ac.version
    gemo.gatts['sensor_lut'] = sensor_lut

    ## call this hDSF
    gemo.gatts['atmospheric_correction_method'] = 'hDSF'
    gemo.gatts['ac_aot_550'] = aer_depth
    for lut in setu['luts']:
        if lut.endswith('MOD1') & (aer_model == 'aerosol_haze 1'):
            aer_lut = lut
        if lut.endswith('MOD2') & (aer_model == 'aerosol_haze 4'):
            aer_lut = lut
        if lut.endswith('MOD3') & (aer_model == 'aerosol_haze 5'):
            aer_lut = lut
    gemo.gatts['ac_model'] = aer_lut

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

    ## compute surface reflectance
    for bi, b in enumerate(gem.bands):
        if ('rhot_ds' not in gem.bands[b]) or ('tt_gas' not in gem.bands[b]):
            if setu['verbosity'] > 2: print('Band {} at {} nm not in bands dataset'.format(b, gem.bands[b]['wave_name']))
            continue

        if gem.bands[b]['rhot_ds'] not in gem.datasets:
            if setu['verbosity'] > 2: print('Band {} at {} nm not in available rhot datasets'.format(b, gem.bands[b]['wave_name']))
            continue ## skip if we don't have rhot for a band that is in the RSR file

        ## input and output  dataset
        dsi = gem.bands[b]['rhot_ds']
        dso = gem.bands[b]['rhos_ds']

        ## exclude/include band if requested
        if setu['l2r_exclude_bands'] is not None:
            if dso in setu['l2r_exclude_bands']:
                print('Skipping {} which is in l2r_exclude_bands'.format(dso))
                continue
        if setu['l2r_include_bands'] is not None:
            if dso not in setu['l2r_include_bands']:
                print('Skipping {} which is not in l2r_include_bands'.format(dso))
                continue

        rhot, cur_att = gem.data(dsi, attributes = True)

        ## store rhot in output file
        if copy_rhot:
           gemo.write(dsi, rhot, ds_att = cur_att)

        ## skip low transmittance bands
        if gem.bands[b]['tt_gas'] < setu['min_tgas_rho']:
            if setu['verbosity'] > 2: print('Band {} at {} nm has tgas < min_tgas_rho ({:.2f} < {:.2f})'.format(b, gem.bands[b]['wave_name'], gem.bands[b]['tt_gas'], setu['min_tgas_rho']))
            continue

        ## remove path reflectance
        rhot_noatm = rhot - lrr[lira_path][b]

        ## correct transmittance (and sa)
        if not lira_sa:
            rhos = (rhot_noatm) / (lrr['Td_tot'][b] * lrr['Tu_tot'][b])
        else:
            rhos = (rhot_noatm) / (lrr['Td_tot'][b] * lrr['Tu_tot'][b] + lrr['sa'][b] * rhot_noatm)

        ## store a/c parameters
        for k in lrr: cur_att[k] = lrr[k][b]

        ## write rhos
        gemo.write(dso, rhos, ds_att = cur_att)

    ## update dataset info
    gemo.setup()

    ## glint correction - should do using libRadtran results
    if setu['dsf_residual_glint_correction']:
        if hyper:
            print('hDSF glint correction not yet implemented for hyperspectral sensors')
        else:
            if setu['dsf_residual_glint_correction_method'] != 'default':
                print('dsf_residual_glint_correction_method={} not implemented after RAdCor'.format(setu['dsf_residual_glint_correction_method']))
            else:
                print('Running glint correction!')
                ## run glint correction
                ret = ac.glint.default(gemo, settings = setu, new_file = False, write = True)

    ## close files
    gem, gemo = None, None

    return(ofile, setu)
