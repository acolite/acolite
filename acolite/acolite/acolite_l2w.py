## def acolite_l2w
## new L2W parameter computation for L2R generic extracted miniscene
## written by Quinten Vanhellemont, RBINS
## 2021-03-09
## modifications: 2021-12-08 (QV) added nc_projection
##                2022-09-28 (QV) changed gem from dict to object
##                2024-01-26 (QV) moved flagging to a separate function
##                2024-03-14 (QV) update settings handling
##                2024-04-16 (QV) changed output NetCDF handling to using gem
##                2025-01-16 (QV) added rrs and lowercase to parameter loop
##                2025-02-04 (QV) added additional S2C_MSI outputs
##                                updated settings handling, added rsr_version, fixed rhos_ds selection
##                2025-03-12 (QV) added wildcards for copy_datasets in separate file
##                2025-03-27 (QV) update to RSR reading

def acolite_l2w(gem, output = None, settings = None,
                target_file = None, load_data = True, new = True, ):

    import os
    import numpy as np
    import acolite as ac
    import scipy.ndimage
    import skimage.color

    ## read gem file if NetCDF
    close = False
    if type(gem) is str:
        gemf = '{}'.format(gem)
        gem = ac.gem.gem(gem)
        close = True
    gemf = gem.file

    ## combine default and user defined settings
    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}
    ## get sensor specific defaults
    setd = ac.acolite.settings.parse(gem.gatts['sensor'])
    ## set sensor default if user has not specified the setting
    for k in setd:
        if k not in ac.settings['user']: setu[k] = setd[k]
    ## end set sensor specific defaults

    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    verbosity = setu['verbosity']

    ## nothing to do if no l2w_parameters requested
    if setu['l2w_parameters'] == None: return(None)

    ## set up output file
    if target_file is None:
        output_name = os.path.basename(gemf).replace('.nc', '')
        output_name = output_name.replace('_L2R', '_L2W')
        if 'output' in setu: output = setu['output']
        odir = output if output is not None else os.path.dirname(gemf)
        ofile = '{}/{}.nc'.format(odir, output_name)
    else:
        ofile = '{}'.format(target_file)

    ## keep data in memory
    gem.store = setu['l2w_data_in_memory']
    if verbosity > 3: print('L2W keeping data in memory: {}'.format(gem.store))
    gem.verbosity = setu['verbosity']

    ## get rhot and rhos wavelengths
    rhot_ds = [ds for ds in gem.datasets if ds.startswith('rhot_')]
    rhot_waves = [int(ds.split('_')[-1]) for ds in rhot_ds]
    if len(rhot_waves) == 0: print('{} is probably not an ACOLITE L2R file: {} rhot datasets.'.format(gemf, len(rhot_ds)))
    if verbosity > 3:
        print('rhot datasets: ', rhot_ds)
        print('rhot wavelengths: ', rhot_waves)

    rhos_ds = [ds for ds in gem.datasets if ds.startswith('rhos_')]
    rhos_waves = [int(ds.split('_')[-1]) for ds in rhos_ds]
    if len(rhos_waves) == 0: print('{} is probably not an ACOLITE L2R file: {} rhos datasets.'.format(gemf, len(rhos_ds)))
    if verbosity > 3:
        print('rhos datasets: ', rhos_ds)
        print('rhos wavelengths: ', rhos_waves)

    if gem.gatts['acolite_file_type'] != 'L2R':
        print('Only L2W processing of ACOLITE L2R files supported.')
        print('{} is a "{}" file'.format(gemf, gem.gatts['acolite_file_type']))
        return(None)

    ## read rsr
    hyper = gem.gatts['sensor'] in ac.config['hyper_sensors']
    lut_sensor = '{}'.format(gem.gatts['sensor'])
    if setu['rsr_version'] is not None: lut_sensor = '{}_{}'.format(gem.gatts['sensor'], setu['rsr_version'])

    ## create or load RSR copied from L2R - maybe to split off in a function
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
        lut_sensor = '{}'.format(gem.gatts['sensor'])
    elif (sensor_version is not None) & (sensor_version in rsrd):
        print('Using RSR {} for sensor {}'.format(sensor_version, gem.gatts['sensor']))
        lut_sensor = '{}'.format(sensor_version)
    else:
        print('Could not find {} RSR'.format(gem.gatts['sensor']))
        return()

    ## spectral turbidity/spm
    nechad_range = setu['nechad_range']
    for k in ['tur_nechad2009_*', 'tur_nechad2009ave_*',
              'spm_nechad2010_*', 'spm_nechad2010ave_*',
              'tur_dogliotti2022_*']:
        if (k in setu['l2w_parameters']):
            setu['l2w_parameters'].remove(k)
            ## add key to auto grouping for SNAP
            grouping_key = k[0:-2]
            grouping_key = grouping_key[0:5].upper() + grouping_key[5:]
            if 'ave' in grouping_key: grouping_key = grouping_key.replace('ave','Ave')
            if grouping_key not in gem.gatts['auto_grouping']:
                gem.gatts['auto_grouping'] += ':{}'.format(grouping_key)
            ## run through bands
            for b in rsrd[lut_sensor]['wave_name']:
                w = rsrd[lut_sensor]['wave_nm'][b]
                wn = rsrd[lut_sensor]['wave_name'][b]
                if (w < nechad_range[0]) or (w > nechad_range[1]): continue
                if wn == 'nan': continue
                setu['l2w_parameters'].append(k.replace('_*', '_{}').format(wn))

    ## determine chl_ocx dataset for chlor_a
    if 'chlor_a' in setu['l2w_parameters']:
        ## remove chlor_a parameter
        setu['l2w_parameters'].remove('chlor_a')
        ## load config
        cx_file = ac.config['data_dir']+'/Shared/algorithms/chl_ci/chl_ocx.txt'
        cx_cfg = ac.acolite.settings.read(cx_file)
        ## determine parameters
        if gem.gatts['sensor'] not in cx_cfg:
            print('chl_ocx not configured for {}'.format(gem.gatts['sensor']))
            print('Not outputting chlor_a')
        else:
            cx_par = cx_cfg[gem.gatts['sensor']]
            print('Using {} as chl_ocx for {}'.format(cx_par, gem.gatts['sensor']))
            if cx_par not in setu['l2w_parameters']:
                setu['l2w_parameters'].append(cx_par)
            if 'chlor_a' not in setu['l2w_parameters']:
                setu['l2w_parameters'].append('chlor_a')
            ## move chl_ci to the end if it is in the list
            setu['l2w_parameters'].append('chl_ci')
            setu['l2w_parameters'].remove('chl_ci')
    ## end determine chl_ocx

    ## compute flag value to mask for water products
    flag_value = 0
    if setu['l2w_mask']:
        flag_value = 2**setu['flag_exponent_swir']
        if setu['l2w_mask_cirrus']: flag_value += 2**setu['flag_exponent_cirrus']
        if setu['l2w_mask_high_toa']: flag_value += 2**setu['flag_exponent_toa']
        if setu['l2w_mask_negative_rhow']: flag_value += 2**setu['flag_exponent_negative']
        if setu['l2w_mask_mixed']: flag_value += 2**setu['flag_exponent_mixed']
        if setu['dem_shadow_mask']: flag_value += 2**setu['flag_exponent_dem_shadow']

    ## compute flags
    l2_flags = ac.acolite.acolite_flags(gem)

    ## read parameters that are copied from L2R -> L2W
    copied_parameters = ac.acolite.settings.read_list(ac.config['data_dir']+'/ACOLITE/l2w_parameters_copy.txt')

    ## set up output file gem
    gemo = ac.gem.gem(ofile, new = new)
    gemo.gatts = {k: gem.gatts[k] for k in gem.gatts}
    gemo.nc_projection = gem.nc_projection
    gemo.gatts['acolite_file_type'] = 'L2W'
    gemo.gatts['ofile'] = ofile
    gemo.verbosity = 6

    ## list datasets to copy over from L2R
    copy_datasets = ['lon', 'lat']
    if verbosity > 3: print('Copying datasets from L2R.')
    for cur_par in gem.datasets:
        if cur_par in copy_datasets: continue

        ## add rhow / Rrs from rhos
        if ('rhos_' in cur_par):
            if (('rhow_*' in setu['l2w_parameters'])) |\
                (cur_par.replace('rhos_', 'rhow_') in setu['l2w_parameters']):
                new_par = cur_par.replace('rhos_', 'rhow_')
                if new_par not in copy_datasets: copy_datasets.append(new_par)
            if (('Rrs_*' in setu['l2w_parameters'])) |\
                (cur_par.replace('rhos_', 'Rrs_') in setu['l2w_parameters']):
                new_par = cur_par.replace('rhos_', 'Rrs_')
                if new_par not in copy_datasets: copy_datasets.append(new_par)
            if (('rrs_*' in setu['l2w_parameters'])) |\
                (cur_par.replace('rhos_', 'rrs_') in setu['l2w_parameters']):
                new_par = cur_par.replace('rhos_', 'rrs_')
                if new_par not in copy_datasets: copy_datasets.append(new_par)

        ## copy parameter if listed in l2w_parameters
        if (cur_par in setu['l2w_parameters']):
            copy_datasets.append(cur_par)
        ## evaluate wildcards, now stored in separate file
        else:
            for copy_par in copied_parameters:
                if (copy_par in setu['l2w_parameters']):
                    base_par = copy_par.replace('*','')
                    if (cur_par.startswith(base_par)) |\
                       (cur_par.lower().startswith(base_par)):
                        copy_datasets.append(cur_par)
    ## end list datasets to copy over from L2R

    ## copy datasets
    for ci, cur_par in enumerate(copy_datasets):
        factor = 1.0
        cur_tag = '{}'.format(cur_par)
        mask = False
        att_add = {}
        ## copy Rrs/rhow
        if 'Rrs_' in cur_par:
            factor = 1.0/np.pi
            cur_tag = cur_par.replace('Rrs_','rhos_')
            mask = True
            att_add = {'algorithm':'Remote sensing reflectance', 'dataset':'rhos'}
            att_add['reference']=''
            att_add['algorithm']=''
        if 'rrs_' in cur_par:
            factor = 1.0/np.pi
            cur_tag = cur_par.replace('rrs_','rhos_')
            mask = True
            att_add = {'algorithm':'Subsurface remote sensing reflectance', 'dataset':'rhos'}
            att_add['reference']=''
            att_add['algorithm']=''
        if 'rhow_' in cur_par:
            factor = 1.0
            cur_tag = cur_par.replace('rhow_','rhos_')
            mask = True
            att_add = {'algorithm':'Water reflectance', 'dataset':'rhos'}
            att_add['reference']=''
            att_add['algorithm']=''

        # skip if the source dataset is not in the inputfile
        if cur_tag not in gem.datasets:
            print('{} not in {}'.format(cur_tag, gemf))
            continue

        ## if data already read copy here
        print('Copying {}, base dataset {}'.format(cur_par, cur_tag))
        cur_data, cur_att = gem.data(cur_tag, attributes=True)
        if factor != 1.0: cur_data *= factor

        ## Lee et al 2002 subsurface reflectance
        if cur_par[0:3] == 'rrs': cur_data /= (0.52 + 1.7 * cur_data)

        ## apply mask to Rrs and rhow
        if (mask) & (setu['l2w_mask_water_parameters']): cur_data[(l2_flags & flag_value)!=0] = np.nan
        if verbosity > 1: print('Writing {}'.format(cur_par))
        ## add attributes
        for k in att_add: cur_att[k] = att_add[k]
        gemo.write(cur_par, cur_data, ds_att = cur_att)

    ## write l2 flags
    gemo.write('l2_flags', l2_flags)

    qaa_computed, p3qaa_computed = False, False

    ## parameter loop
    ## compute other parameters
    for cur_par in setu['l2w_parameters']:
        if cur_par in copied_parameters: continue
        if cur_par.lower() in copied_parameters + ['rrs_*', 'rhow_*', '', ' ']: continue ## we have copied these above
        if cur_par.lower() in [ds.lower() for ds in ac.shared.nc_datasets(ofile)]: continue ## parameter already in output dataset (would not work if we are appending subsets to the ncdf)
        if cur_par.lower().startswith('bt'): continue
        cur_par = cur_par.lower() ## workaround for parameters provided with upper case

        ## split on underscores
        sp = cur_par.split('_')

        ## default mask and empty dicts for current parameter
        mask = False
        par_data = {}
        par_atts = {}

        #############################
        ## Nechad turbidity/spm
        if 'nechad' in cur_par:
            mask = True ## water parameter so apply mask

            nechad_parameter = 'centre'
            if '2016' in cur_par: nechad_parameter = '2016'
            elif 'ave' in cur_par: nechad_parameter = 'resampled'

            ## turbidity
            if cur_par[0] == 't':
                nechad_par = 'TUR'
                par_attributes = {'algorithm':'Nechad et al. 2009', 'title':'Nechad Turbidity'}
                par_attributes['reference']='Nechad et al. 2009'
                par_attributes['algorithm']='2009 calibration'
                cal_year = 2009
            elif cur_par[0] == 's':
                nechad_par = 'SPM'
                par_attributes = {'algorithm':'Nechad et al. 2010', 'title':'Nechad SPM'}
                par_attributes['reference']='Nechad et al. 2010'
                par_attributes['algorithm']='2010 calibration'
                cal_year = 2010
            else:
                continue

            ## find out which band to use
            nechad_band = None
            nechad_wave = 665
            sp = cur_par.split('_')
            if 'ave' not in cur_par:
                if len(sp) > 2: nechad_wave = sp[2]
            else:
                if len(sp) > 2: nechad_wave = sp[2]
                if len(sp) > 3: nechad_wave = sp[3]

            ## find out wavelength for turbidity product
            if type(nechad_wave) == str:
                if nechad_wave.lower() == 'red': nechad_wave=665
                elif nechad_wave.lower() == 'nir': nechad_wave=865
                elif nechad_wave.lower() == 'green': nechad_wave=560
                else:
                    try: ## if wavelength is specified
                        nechad_wave = int(nechad_wave)
                    except:
                        continue

            ci, cw = ac.shared.closest_idx(rhos_waves, int(nechad_wave))
            for b in rsrd[lut_sensor]['rsr_bands']:
                if (rsrd[lut_sensor]['wave_name'][b] == str(cw)): nechad_band = b

            A_Nechad, C_Nechad = None, None

            ## band center
            if nechad_parameter == 'centre':
                par_base = '{}_Nechad{}'.format(nechad_par, cal_year)
                A_Nechad, C_Nechad = ac.parameters.nechad.coef_band(nechad_par, cal_year = cal_year, wave = cw)

            ## resampled to band
            if nechad_parameter == 'resampled':
                par_base = '{}_Nechad{}Ave'.format(nechad_par, cal_year)
                A_Nechad, C_Nechad = ac.parameters.nechad.coef_band(nechad_par, cal_year = cal_year, \
                                                                sensor = gem.gatts['sensor'], band = nechad_band, \
                                                                rsr = rsrd[lut_sensor]['rsr'])

            ## resampled to band by Nechad 2016
            if nechad_parameter == '2016':
                par_base = '{}_Nechad2016'.format(nechad_par)
                par_attributes['algorithm']='2016 calibration'
                nechad_dict = ac.parameters.nechad.coef_2016()
                for k in nechad_dict:
                    if k['sensor'] != gem.gatts['sensor'].split('_')[-1]: continue
                    if k['band'] != 'B{}'.format(nechad_band): continue
                    if k['par'] != nechad_par: continue
                    A_Nechad = k['A']
                    C_Nechad = k['C']

            ## if we have A and C we can continue
            if (A_Nechad is not None) & (C_Nechad is not None):
                #cur_ds = 'rhos_{}'.format(cw)
                cur_ds = [ds for ds in rhos_ds if ('{:.0f}'.format(cw) in ds)][0]
                par_name = '{}_{}'.format(par_base, cur_ds.replace('rhos_',''))
                cur_data = 1.0 * gem.data(cur_ds)
                ## compute parameter
                cur_mask = np.where(cur_data >= (setu['nechad_max_rhow_C_factor'] * C_Nechad))
                cur_data = (A_Nechad * cur_data) / (1.-(cur_data/C_Nechad))
                cur_data[cur_mask] = np.nan
                par_data[par_name] = cur_data
                par_atts[par_name] = par_attributes
                par_atts[par_name]['ds_name'] = par_name
                par_atts[par_name]['A_{}'.format(nechad_par)] = A_Nechad
                par_atts[par_name]['C_{}'.format(nechad_par)] = C_Nechad
        ## end nechad turbidity/spm
        #############################

        #############################
        ## start Dogliotti 2015 turbidity
        if 'dogliotti2015' in cur_par:
            mask = True ## water parameter so apply mask

            par_attributes = {'algorithm':'Dogliotti et al. 2015', 'title':'Dogliotti Turbidity'}
            par_attributes['reference']='Dogliotti et al. 2015'
            par_attributes['algorithm']=''

            ## get config
            par_name = 'TUR_Dogliotti2015'
            dcfg = 'defaults'
            dogliotti_par = 'blended'
            if len(sp) > 2:
                if sp[2] not in ['red', 'nir']:
                    dcfg = sp[2]
                else:
                    dogliotti_par = sp[2]
                par_name+='_{}'.format(dogliotti_par)
            par_attributes['ds_name'] = par_name
            cfg = ac.parameters.dogliotti.coef(config=dcfg)

            ## identify bands
            ri, rw = ac.shared.closest_idx(rhos_waves, int(cfg['algo_wave_red']))
            ni, nw = ac.shared.closest_idx(rhos_waves, int(cfg['algo_wave_nir']))
            for b in rsrd[lut_sensor]['rsr_bands']:
                if (rsrd[lut_sensor]['wave_name'][b] == str(rw)): red_band = b
                if (rsrd[lut_sensor]['wave_name'][b] == str(nw)): nir_band = b

            ## store settings in atts
            for k in cfg: par_attributes[k] = cfg[k]
            par_attributes['wave_red'] = rw
            par_attributes['wave_nir'] = nw

            ## read red data
            #cur_ds = 'rhos_{}'.format(rw)
            cur_ds = [ds for ds in rhos_ds if ('{:.0f}'.format(rw) in ds)][0]
            red = 1.0 * gem.data(cur_ds)
            tur = (par_attributes['A_T_red'] * red) / (1.-red/par_attributes['C_T_red'])

            ## read nir data
            #cur_ds = 'rhos_{}'.format(nw)
            cur_ds = [ds for ds in rhos_ds if ('{:.0f}'.format(nw) in ds)][0]
            nir = 1.0 * gem.data(cur_ds)
            nir_tur = (par_attributes['A_T_nir'] * nir) / (1.-nir/par_attributes['C_T_nir'])

            if dogliotti_par == 'blended':
                ## replace most turbid with nir band
                dsub = np.where(red >= par_attributes['upper_lim'])
                if len(dsub[0]) > 0: tur[dsub] = nir_tur[dsub]
                ## blend in between
                dsub = np.where((red < par_attributes['upper_lim']) & (red >= par_attributes['lower_lim']))
                if len(dsub[0]) > 0:
                    w=(red[dsub]  - par_attributes['lower_lim']) / (par_attributes['upper_lim']-par_attributes['lower_lim'])
                    tur[dsub] = ((1.-w) * tur[dsub]) + (w*nir_tur[dsub])
                par_data[par_name] = tur
                par_atts[par_name] = par_attributes
            elif dogliotti_par == 'red':
                par_data[par_name] = tur
                par_atts[par_name] = par_attributes
            elif dogliotti_par == 'nir':
                par_data[par_name] = nir_tur
                par_atts[par_name] = par_attributes
            red = None
            nir = None
        ## end Dogliotti turbidity
        #############################

        #############################
        ## Dogliotti generic turbidity
        if 'dogliotti2022' in cur_par:
            mask = True ## water parameter so apply mask
            nechad_par = 'TUR'
            par_attributes = {'algorithm':'Dogliotti et al. in prep', 'title':'Dogliotti Turbidity'}
            par_attributes['reference']='Nechad et al. 2009'
            par_attributes['algorithm']='turbidity recalibration 2022'
            nechad_parameter = 'center'
            if not hyper: nechad_parameter = 'resampled'

            ## find out which band to use
            nechad_band = None
            nechad_wave = 665
            sp = cur_par.split('_')
            if len(sp) > 2: nechad_wave = sp[2]
            if len(sp) > 3: nechad_wave = sp[3]

            ## find out wavelength for turbidity product
            if type(nechad_wave) == str:
                if nechad_wave.lower() == 'red': nechad_wave=665
                elif nechad_wave.lower() == 'nir': nechad_wave=865
                elif nechad_wave.lower() == 'green': nechad_wave=560
                else:
                    try: ## if wavelength is specified
                        nechad_wave = int(nechad_wave)
                    except:
                        continue

            ci, cw = ac.shared.closest_idx(rhos_waves, int(nechad_wave))
            for b in rsrd[lut_sensor]['rsr_bands']:
                if (rsrd[lut_sensor]['wave_name'][b] == str(cw)): nechad_band = b

            A_Nechad, C_Nechad = None, None
            nechad_dict = ac.parameters.dogliotti.coef_hyper()

            ## band center
            if nechad_parameter == 'center':
                didx,algwave = ac.shared.closest_idx(nechad_dict['wave'], cw)
                A_Nechad = nechad_dict['A_mean'][didx]
                C_Nechad = nechad_dict['C'][didx]
            elif nechad_parameter == 'resampled':
                # resample parameters to band
                cdct = ac.shared.rsr_convolute_dict(nechad_dict['wave']/1000, nechad_dict['C'], rsrd[lut_sensor]['rsr'])
                adct = ac.shared.rsr_convolute_dict(nechad_dict['wave']/1000, 1/nechad_dict['A_mean'], rsrd[lut_sensor]['rsr'])
                adct = {k:1/adct[k] for k in adct}
                A_Nechad = adct[nechad_band]
                C_Nechad = cdct[nechad_band]

            ## if we have A and C we can continue
            if (A_Nechad is not None) & (C_Nechad is not None):
                cur_ds = [ds for ds in rhos_ds if ('{:.0f}'.format(cw) in ds)][0]
                par_name = '{}_Dogliotti2022_{}'.format(nechad_par, cur_ds.replace('rhos_',''))
                cur_data = 1.0 * gem.data(cur_ds)
                ## compute parameter
                cur_mask = np.where(cur_data >= (setu['nechad_max_rhow_C_factor'] * C_Nechad))
                cur_data = (A_Nechad * cur_data) / (1.-(cur_data/C_Nechad))
                cur_data[cur_mask] = np.nan
                par_data[par_name] = cur_data
                par_atts[par_name] = par_attributes
                par_atts[par_name]['ds_name'] = par_name
                par_atts[par_name]['wavelength'] = cw
                par_atts[par_name]['A_{}'.format(nechad_par)] = A_Nechad
                par_atts[par_name]['C_{}'.format(nechad_par)] = C_Nechad
        ## end Dogliotti generic turbidity
        #############################

        #############################
        ## start Novoa blended turbidity
        if 'novoa2017' in cur_par:
            mask = True ## water parameter so apply mask
            novoa_par = {'t':'TUR', 's':'SPM'}[cur_par.split('_')[0][0].lower()]
            par_attributes = {'algorithm':'Novoa et al. 2017', 'title':'Novoa Turbidity'}
            par_attributes['reference']='Novoa et al. 2017'

            ## read config dict
            novoa_config='defaults'
            novoa_dict = ac.shared.import_config('{}/Shared/algorithms/Novoa/{}.txt'.format(ac.config['data_dir'], novoa_config), parse=True)

            ## find suitable bands for Novoa red-NIR switching algorithm
            tur_waves = []
            tur_bands = []
            for b in rsrd[lut_sensor]['wave_name']:
                w = rsrd[lut_sensor]['wave_nm'][b]
                if (w >= setu['nechad_range'][0]) and (w <= setu['nechad_range'][1]):
                    tur_waves.append(w)
                    tur_bands.append(b)

            ## number of products possible (for VIIRS use both I and M bands)
            novoa_products = 1
            if 'VIIRS' in gem.gatts['sensor']: novoa_products = 2

            ## find red and NIR bands to use
            novoa_bands = []
            novoa_waves = []
            for pi in range(novoa_products):
                novoa_bands.append([])
                novoa_waves.append([])

                for ci, cw in enumerate(novoa_dict['novoa_waves_req']):
                    ## compute Novoa SPM for both I and M bands
                    if 'VIIRS' in gem.gatts['sensor']:
                        band_prefix = {0:'I', 1:'M'}[pi]
                        tur_waves_ = [w for ii, w in enumerate(tur_waves) if tur_bands[ii][0] == band_prefix]
                        tur_bands_ = [b for b in tur_bands if b[0] == band_prefix]
                        wi, wv = ac.shared.closest_idx(tur_waves_, cw)
                        b = tur_bands_[wi]
                    else:
                        wi, wv = ac.shared.closest_idx(tur_waves, cw)
                        b = tur_bands[wi]
                    if wv > [cw+novoa_dict['novoa_waves_off'][ci]]: continue
                    if wv < [cw-novoa_dict['novoa_waves_off'][ci]]: continue
                    novoa_bands[pi].append(b)
                    novoa_waves[pi].append(wv)

            if 'nechad' in novoa_dict['novoa_algorithm']:
                print(novoa_dict)

                ## run through needed products
                for pi in range(novoa_products):
                    if len(novoa_bands[pi]) != 2: continue
                    ## used bands
                    redb, nirb = novoa_bands[pi]
                    redw, nirw = novoa_waves[pi]
                    ## output parameter name
                    par_name = '{}_Novoa2017'.format(novoa_par)
                    if 'VIIRS' in gem.gatts['sensor']:
                        par_name += '_{}'.format(redb[0].upper())

                    ## find proper rhos datasets
                    red_ds = [ds for ds in rhos_ds if ('{:.0f}'.format(redw) in ds)][0]
                    nir_ds = [ds for ds in rhos_ds if ('{:.0f}'.format(nirw) in ds)][0]

                    ## find Nechad calibration for this band
                    A_Nechad_red, C_Nechad_red = None, None
                    A_Nechad_nir, C_Nechad_nir = None, None
                    if novoa_dict['novoa_algorithm'] == 'nechad_centre':
                        A_Nechad_red, C_Nechad_red = ac.parameters.nechad.coef_band(cur_par.split('_')[0], wave=redw)
                        A_Nechad_nir, C_Nechad_nir = ac.parameters.nechad.coef_band(cur_par.split('_')[0], wave=nirw)
                    if novoa_dict['novoa_algorithm'] == 'nechad_average':
                        A_Nechad_red, C_Nechad_red = ac.parameters.nechad.coef_band(cur_par.split('_')[0], sensor=gem.gatts['sensor'], band=redb)
                        A_Nechad_nir, C_Nechad_nir = ac.parameters.nechad.coef_band(cur_par.split('_')[0], sensor=gem.gatts['sensor'], band=nirb)

                    ## if we have A and C we can continue
                    if (A_Nechad_red is not None) & (C_Nechad_red is not None) &\
                       (A_Nechad_nir is not None) & (C_Nechad_nir is not None):
                        par_attributes['A_Nechad_red']=A_Nechad_red
                        par_attributes['C_Nechad_red']=C_Nechad_red
                        par_attributes['A_Nechad_nir']=A_Nechad_nir
                        par_attributes['C_Nechad_nir']=C_Nechad_nir

                        ## load datasets
                        red_data = 1.0 * gem.data(red_ds)
                        nir_data = 1.0 * gem.data(nir_ds)

                        ## compute red and NIR turbidity
                        red_tur = (A_Nechad_red * red_data) / (1.-(red_data/C_Nechad_red))
                        red_tur[red_tur<0] = np.nan
                        nir_tur = (A_Nechad_nir * nir_data) / (1.-(nir_data/C_Nechad_nir))
                        nir_tur[nir_tur<0] = np.nan

                        ## mask Nechad product in non-linear regime
                        #red_mask = np.where(red_data >= (setu['nechad_max_rhow_C_factor'] * C_Nechad_red))
                        #red_tur[red_mask] = np.nan
                        #nir_mask = np.where(nir_data >= (setu['nechad_max_rhow_C_factor'] * C_Nechad_nir))
                        #nir_tur[nir_mask] = np.nan
                        nir_data = None

                        ## determine switching bounds to use red/NIR bands
                        redi = np.where(red_data <= novoa_dict['rhow_switch_nir'][0])
                        niri = np.where(red_data >= novoa_dict['rhow_switch_nir'][1])

                        ## log blending bounds ratio
                        logb = np.log(novoa_dict['rhow_switch_nir'][1]/novoa_dict['rhow_switch_nir'][0])

                        ## blended product
                        tur_blend = red_tur * np.log(novoa_dict['rhow_switch_nir'][1]/red_data) / logb + \
                                    nir_tur * np.log(red_data/novoa_dict['rhow_switch_nir'][0]) / logb

                        ## indices for red-NIR blend
                        rnblend = np.where((red_data > novoa_dict['rhow_switch_nir'][0]) &\
                                           (red_data < novoa_dict['rhow_switch_nir'][1]) &\
                                           (np.isfinite(tur_blend)))

                        ## composite algorithm
                        par_data[par_name] = red_tur
                        par_data[par_name][niri] = nir_tur[niri]
                        par_data[par_name][rnblend] = tur_blend[rnblend]
                        par_atts[par_name] = par_attributes

                        ## track algorithm
                        if novoa_dict['novoa_output_switch']:
                            alg = np.zeros(par_data[par_name].shape)
                            alg[redi] = 3
                            alg[rnblend] = 4
                            alg[niri] = 5
                            par_data[par_name+'_switch'] = alg
                            par_atts[par_name+'_switch'] = par_attributes
            else:
                print('Algorithm "{}" for Novoa not supported.'.format(novoa_dict['novoa_algorithm']))
                continue
        ## end Novoa blended turbidity
        #############################

        #############################
        ## CHL_OC
        if 'chl_oc' in cur_par:
            mask = True ## water parameter so apply mask
            ## load config
            chl_oc_wl_diff = 20
            cfg = ac.parameters.chl_oc.coef()
            if gem.gatts['sensor'] not in cfg:
                print('{} not configured for {}'.format(cur_par, gem.gatts['sensor']))
                continue

            par_attributes = {'algorithm':'Chlorophyll a blue/green ratio', 'dataset':'rhos'}
            par_attributes['reference']='Franz et al. 2015'
            par_attributes['algorithm']=''

            if (cur_par == 'chl_oc') | (cur_par == 'chl_oc2'):
                par_name = 'chl_oc2'
            elif (cur_par == 'chl_oc3'):
                par_name = 'chl_oc3'
            elif (cur_par == 'chl_oc4'):
                par_name = 'chl_oc4'
            else:
                par_name = cur_par

            if par_name not in cfg[gem.gatts['sensor']]:
                print('{} not configured for {}'.format(par_name, gem.gatts['sensor']))
                continue
            par_attributes['ds_name']=par_name
            chl_dct = cfg[gem.gatts['sensor']][par_name]

            ## get bands
            blue, green = None, None
            for w in chl_dct['blue'] + chl_dct['green']:
                ci, cw = ac.shared.closest_idx(rhos_waves, int(w))
                if np.abs(cw-w) > chl_oc_wl_diff: continue
                cur_data = 1.0 * gem.data(rhos_ds[ci])
                if w in chl_dct['blue']:
                    if blue is None:
                        par_attributes['blue_wave_sel'] = [cw]
                        blue = 1.0 * cur_data
                    else:
                        par_attributes['blue_wave_sel'] += [cw]
                        blue = np.nanmax((blue, cur_data), axis=0)
                if w in chl_dct['green']:
                    if green is None:
                        par_attributes['green_wave_sel'] = [cw]
                        green = 1.0 * cur_data
            if (blue is None) | (green is None): continue

            ## compute CHL
            ratio = np.log10(blue/green)
            blue, green = None, None
            par_data[par_name] = chl_dct['chl_coef'][0] + chl_dct['chl_coef'][1] * ratio + \
                                                     chl_dct['chl_coef'][2] * ratio * ratio + \
                                                     chl_dct['chl_coef'][3] * ratio * ratio * ratio + \
                                                     chl_dct['chl_coef'][4] * ratio * ratio * ratio * ratio
            par_data[par_name] = np.power(10, par_data[par_name])
            par_atts[par_name] = par_attributes
        ## end CHL_OC
        #############################

        #############################
        ## CHL_CI
        if (cur_par[0:6] == 'chl_ci') | (cur_par[0:7] == 'chlor_a'):
            mask = True ## water parameter so apply mask
            ## load config
            ci_file = ac.config['data_dir']+'/Shared/algorithms/chl_ci/defaults.txt'
            ci_cfg = ac.acolite.settings.read(ci_file)
            for k in ci_cfg:
                if type(ci_cfg[k]) == bool: continue
                if type(ci_cfg[k]) == list:
                    ci_cfg[k] = [float(v) for v in ci_cfg[k]]
                else:
                    ci_cfg[k] = float(ci_cfg[k])
            par_attributes = {'algorithm':'Chlorophyll a CI', 'dataset':'rhos'}
            par_attributes['reference']='Hu et al. 2012'

            ## find required datasets
            required_datasets = []
            selected_waves = []
            green_diff = 0
            for w in ci_cfg['wave']:
                wi, ws = ac.shared.closest_idx(rhos_waves, w)
                if w == ci_cfg['wave'][1]:
                    green_wave = ws
                    green_diff = np.abs(ws-w)
                required_datasets.append(rhos_ds[wi])
                selected_waves.append(float(ws))

            ## determine green shift
            green_shift = None
            if green_diff > ci_cfg['wave_green_diff']:
                for k in ci_cfg:
                    if 'shift' in k:
                        if (green_wave >= ci_cfg[k][0]) & (green_wave <= ci_cfg[k][1]):
                            green_shift = ci_cfg[k]
                            selected_waves[1] = ci_cfg['wave'][1]
                            break
                if green_shift == None:
                    print('Green band cannot be shifted for CI')
                    continue

            ## get data
            for di, cur_ds in enumerate(required_datasets):
                print(cur_ds)
                if di == 0: tmp_data = []
                cur_data = gem.data(cur_ds) / np.pi
                tmp_data.append(cur_data)

            ## perform green shift
            if green_shift != None:
                ## determine pixels for both shifting methods
                sub1 = np.where(tmp_data[1] < green_shift[2])
                sub2 = np.where(tmp_data[1] >= green_shift[2])
                ## shift data
                tmp_data[1][sub1] = np.power(10, (green_shift[3]*np.log10(tmp_data[1][sub1])-green_shift[4]))
                tmp_data[1][sub2] = green_shift[5] * tmp_data[1][sub2] + green_shift[6]

            ## compute ci
            par_data['ci'] = tmp_data[1] - (tmp_data[0] \
                                + (selected_waves[1]-selected_waves[0]) / (selected_waves[2] - selected_waves[0]) \
                                * (tmp_data[2]-tmp_data[0]))
            par_atts['ci'] = par_attributes
            tmp_data = None

            ## compute chl_ci
            par_data['chl_ci'] = np.power(10, ci_cfg['a0CI'] + ci_cfg['a1CI'] * par_data['ci'])
            par_atts['chl_ci'] = par_attributes

            ## compute chlor_a
            if (cur_par[0:7] == 'chlor_a') | ('chlor_a' in setu['l2w_parameters']):
                ## read chl_ocx
                par_data['chl_ocx'], par_atts['chl_ocx'] = ac.shared.nc_data(ofile, cx_par, attributes=True)

                ## copy chl_ci
                par_data['chlor_a'] = par_data['chl_ci'] * 1.0
                par_atts['chlor_a'] = par_attributes
                ## blend chl_ocx
                sub = (par_data['chl_ci'] > ci_cfg['t1']) & (par_data['chl_ci'] <= ci_cfg['t2'])
                alpha = (par_data['chl_ci'][sub] - ci_cfg['t1']) / (ci_cfg['t2']-ci_cfg['t1'])
                beta = (ci_cfg['t2'] - par_data['chl_ci'][sub]) / (ci_cfg['t2']-ci_cfg['t1'])
                par_data['chlor_a'][sub] = alpha*par_data['chl_ocx'][sub] + beta*par_data['chl_ci'][sub]
                ## replace chl_ocx
                sub = par_data['chl_ci'] > ci_cfg['t2']
                par_data['chlor_a'][sub] = par_data['chl_ocx'][sub]

            ## mask chl_ci outside of bounds
            if ci_cfg['mask_t2']:
                sub = par_data['chl_ci'] > ci_cfg['t2']
                par_data['chl_ci'][sub] = np.nan
        ## end CHL_CI
        #############################

        #############################
        ## CHL_RE
        if (cur_par[0:6] == 'chl_re'):
            par_name = cur_par
            mask = True ## apply non water mask
            if gem.gatts['sensor'] not in ['S2A_MSI', 'S2B_MSI', 'S2C_MSI', 'S3A_OLCI', 'S3B_OLCI', 'EN1_MERIS'] + ac.config['hyper_sensors']:
                print('Parameter {} not configured for {}.'.format(par_name,gem.gatts['sensor']))
                continue

            par_split = par_name.split('_')
            par_attributes = {'algorithm':'Red-edge chlorophyll', 'dataset':'rhos'}
            par_attributes['reference']=''
            par_attributes['algorithm']=''

            ## find algorithms
            if len(par_split) >= 3:
                chl_re_algorithm = None

                ###########
                ## GONS
                if par_split[2][0:4] == 'gons':
                    chl_re_algorithm = 'gons'
                    par_attributes['algorithm']='Gons et al. 3 band'
                    par_attributes['reference']='Gons et al. 2005'
                    gons = ac.parameters.chl_re.coef_gons()

                    if par_split[2] == 'gons':
                        req_waves = [670,705,780]
                        gons_name = 'chl_re_gons'
                    elif par_split[2] == 'gons740':
                        req_waves = [670,705,740]
                        gons_name = 'chl_re_gons740'
                    else:
                        print('Parameter {} not configured for {}.'.format(par_name,gem.gatts['sensor']))
                        continue
                    req_waves = [gons[gons_name][tag] for tag in ['red_band', 'rededge_band', 'nir_band']]
                ## end GONS
                ###########

                ###########
                ## MOSES
                if par_split[2][0:5] == 'moses':
                    chl_re_algorithm = 'moses'
                    par_attributes['reference']='Moses et al. 2012'
                    par_attributes['algorithm']='Moses et al. 3 band'
                    par_attributes['a']=(232.29,23.173) ## put coefficients in external file
                    ### get required datasets
                    if par_split[2] == 'moses3b':
                        req_waves = [670,705,780]
                    elif par_split[2] == 'moses3b740':
                        req_waves = [670,705,740]
                    else:
                        print('Parameter {} not configured for {}.'.format(par_name,gem.gatts['sensor']))
                        continue
                ## end MOSES
                ###########

                ###########
                ## MISHRA
                if par_split[2][0:6] == 'mishra':
                    chl_re_algorithm = 'mishra'
                    par_attributes['reference']='Mishra et al. 2014'
                    par_attributes['algorithm']='Mishra et al. 2014, NDCI'
                    par_attributes['a']=(14.039, 86.115, 194.325) ## put coefficients in external file
                    ### get required datasets
                    req_waves = [670,705]
                ## end MISHRA
                ###########

                ###########
                ## BRAMICH
                if par_split[2][0:7] == 'bramich':
                    chl_re_algorithm = 'bramich'
                    par_attributes['reference']='Bramich et al. 2021'
                    par_attributes['algorithm']='Bramich et al. 2021'
                    ### get required datasets
                    req_waves = [670,705,780]
                ## end BRAMICH
                ###########

                ###########
                ## read data
                if chl_re_algorithm is not None:
                    required_datasets,req_waves_selected = [],[]
                    ds_waves = [w for w in rhos_waves]
                    for i, reqw in enumerate(req_waves):
                        widx,selwave = ac.shared.closest_idx(ds_waves, reqw)
                        if abs(float(selwave)-float(reqw)) > 10: continue
                        required_datasets.append(rhos_ds[widx])
                        req_waves_selected.append(selwave)
                    par_attributes['waves']=req_waves_selected
                    if len(required_datasets) != len(req_waves):
                        print('Required datasets not found for {}.'.format(par_name))
                        continue
                    ## get data
                    for di, cur_ds in enumerate(required_datasets):
                        if di == 0: tmp_data = []
                        cur_data = 1.0 * gem.data(cur_ds)
                        tmp_data.append(cur_data)
                else:
                    print('Parameter {} not configured for {}.'.format(par_name,gem.gatts['sensor']))
                    continue
                ## end read data
                ###########

                ###########
                ## compute GONS
                if chl_re_algorithm == 'gons':
                    for k in gons[gons_name]: par_attributes[k] = gons[gons_name][k]
                    gc = gons[gons_name]['chl_coef']
                    bb = (gc[0] * tmp_data[2]) / (gc[1] - gc[2] * tmp_data[2])
                    rm = tmp_data[1]/tmp_data[0]
                    par_data[par_name] = ((rm * (gc[3] + bb)) - gc[4] - np.power(bb,gc[5]))
                    par_data[par_name] /= gons[gons_name]['astar_chl']
                    if '_nocheck' not in par_name:
                        par_data[par_name][par_data[par_name]<0]=np.nan
                        gm = gons[gons_name]['validity']
                        par_data[par_name][((tmp_data[0] <= gm[0]) | (tmp_data[1]/tmp_data[0] <= gm[1]))]=np.nan
                    par_atts[par_name] = par_attributes
                    tmp_data = None
                ## end compute GONS
                ###########

                ###########
                ## compute MOSES
                if chl_re_algorithm == 'moses':
                    par_data[par_name] = par_attributes['a'][0]* \
                                        ((np.power(tmp_data[0],-1)-np.power(tmp_data[1],-1)) * \
                                        tmp_data[2]) + par_attributes['a'][1]
                    if '_nocheck' not in par_name:
                        par_data[par_name][par_data[par_name]<0]=np.nan
                    par_atts[par_name] = par_attributes
                    tmp_data = None
                ## end compute MOSES
                ###########

                ###########
                ## compute MISHRA
                if chl_re_algorithm == 'mishra':
                    ndci = (tmp_data[1]-tmp_data[0]) / (tmp_data[1]+tmp_data[0])
                    par_data[par_name] = par_attributes['a'][0] + par_attributes['a'][1]*ndci + par_attributes['a'][2]*ndci*ndci
                    ndci = None
                    if '_nocheck' not in par_name:
                        par_data[par_name][par_data[par_name]<0]=np.nan
                    par_atts[par_name] = par_attributes
                    tmp_data = None
                ## end compute MISHRA
                ###########

                ###########
                ## compute BRAMICH
                if chl_re_algorithm == 'bramich':
                    ## same as gons
                    bb = (1.61 * tmp_data[2]) / (0.082 - 0.6 * tmp_data[2])
                    rm = tmp_data[1]/tmp_data[0]
                    par_data[par_name] = ((rm * (0.7 + bb)) - 0.40 - np.power(bb,1.05))
                    ## from eq 12 in Bramich 2021
                    #par_data[par_name] = np.power(par_data[par_name]/0.022,1.201)
                    ## from Bramich via RG
                    par_data[par_name] = np.power(par_data[par_name],1.1675)/0.0109
                    if '_nocheck' not in par_name:
                        par_data[par_name][par_data[par_name]<0]=np.nan
                    par_atts[par_name] = par_attributes
                    tmp_data = None
                ## end compute BRAMICH
                ###########
        ## end CHL_RE
        #############################

        #############################
        ## QAA
        if (cur_par[0:3] == 'qaa') | (cur_par == 'qaa5') | (cur_par == 'qaa6') | (cur_par == 'qaaw') |\
           ((cur_par[0:3] == 'qaa') & (('_v5' in cur_par) | ('_v6' in cur_par) | ('_vw' in cur_par))):
            if qaa_computed: continue
            print('QAA')
            mask = True ## water parameter so apply mask
            sensor = gem.gatts['sensor']
            if sensor not in ['L8_OLI', 'L9_OLI', 'S2A_MSI', 'S2B_MSI', 'S2C_MSI']:
                print('QAA not configured for {}'.format(gem.gatts['sensor']))
                continue

            par_attributes = {'algorithm':'QAA', 'dataset':'rhos'}
            par_attributes['standard_name']='qaa'
            par_attributes['long_name']='Quasi Analytical Algorithm outputs'
            par_attributes['units']='various'
            par_attributes['reference']='Lee et al. 2002'
            par_attributes['algorithm']=''

            qaa_coef = ac.parameters.qaa.qaa_coef()
            qaa_wave = [443, 490, 560, 665]
            sen_wave = []

            for ki, k in enumerate(qaa_wave):
                ci, cw = ac.shared.closest_idx(rhos_waves, k)
                #cur_ds = 'rhos_{}'.format(cw)
                cur_ds = [ds for ds in rhos_ds if ('{:.0f}'.format(cw) in ds)][0]
                sen_wave.append(cw)
                if ki == 0: qaa_in = {}
                cur_data = 1.0 * gem.data(cur_ds)
                ## mask data
                if (mask) & (setu['l2w_mask_water_parameters']): cur_data[(l2_flags & flag_value)!=0] = np.nan
                ## convert to Rrs
                qaa_in[k] = cur_data/np.pi

            ## get sun zenith angle
            if 'sza' in gem.datasets:
                sza = gem.data('sza')
            else:
                sza = gem.gatts['sza']

            ## run qaa
            ret = ac.parameters.qaa.qaa_compute(qaa_in, qaa_coef = qaa_coef,
                                        sza=sza, satellite=sensor[0:2])
            ## list possible output parameters
            qaa_pars = list(ret.keys())
            cur_par_out = []
            ## check which parameters are wanted
            cur_par_out = []
            for p in setu['l2w_parameters']:
                for qp in qaa_pars:
                    k = 'qaa_{}'.format(qp)
                    if k in cur_par_out: continue
                    if p == 'qaa':
                        cur_par_out.append(k)
                    elif p in ['qaa5', 'qaa6', 'qaaw']:
                        if (('v5_' not in k) & ('v6_' not in k) & ('vw_' not in k)) |\
                            ('v{}_'.format(p[-1]) in k):
                            cur_par_out.append(k)
                    else:
                        if p == k: cur_par_out.append(k)
            cur_par_out.sort()

            ## reformat for output
            for p in cur_par_out:
                par_data[p] = ret[p[4:]] * 1.0
                ret[p[4:]] = None
                par_atts[p] = par_attributes
                par_atts[p]['ds_name'] = p
            ret = None
            qaa_computed = True
        ## end QAA
        #############################

        #############################
        ## Pitarch 3 band QAA
        if cur_par[0:5] == 'p3qaa':
            if p3qaa_computed: continue

            mask = True ## water parameter so apply mask
            ## load config
            cfg = ac.parameters.pitarch.p3qaa_coef()
            if gem.gatts['sensor'] not in cfg:
                print('P3QAA not configured for {}'.format(gem.gatts['sensor']))
                continue

            par_attributes = {'algorithm':'Pitarch and Vanhellemont, submitted'}

            ## read Blue Green Red data, convert to Rrs
            for k in ['B', 'G', 'R']:
                ci, cw = ac.shared.closest_idx(rhos_waves, cfg[gem.gatts['sensor']]['center_wl'][k])
                #cur_ds = 'rhos_{}'.format(cw)
                cur_ds = [ds for ds in rhos_ds if ('{:.0f}'.format(cw) in ds)][0]

                cur_data = 1.0 * gem.data(cur_ds)
                if (mask) & (setu['l2w_mask_water_parameters']): cur_data[(l2_flags & flag_value)!=0] = np.nan
                if k == 'B':
                    B = cur_data / np.pi
                    Bwave = cw
                if k == 'G':
                    G = cur_data / np.pi
                    Gwave = cw
                if k == 'R':
                    R = cur_data / np.pi
                    Rwave = cw
                cur_data = None

            ## compute 3 band QAA
            print('Computing Pitarch 3 band QAA')
            print('R{} G{} B{}'.format(Rwave, Gwave, Bwave))
            ret = ac.parameters.pitarch.p3qaa_compute(gem.gatts['sensor'], B, G, R)

            ## list possible output parameters
            p3_pars = []
            for k in ret:
                if len(ret[k].shape) == 3:
                    p3_pars += ['p3qaa_{}_{}'.format(k, Bwave), 'p3qaa_{}_{}'.format(k, Gwave), 'p3qaa_{}_{}'.format(k, Rwave)]
                else:
                    p3_pars += ['p3qaa_{}'.format(k)]
            ## check which parameters are wanted
            if 'p3qaa' in setu['l2w_parameters']:
                cur_par_out = p3_pars
            else:
                cur_par_out = [k for k in setu['l2w_parameters'] if k in p3_pars]
                for k in setu['l2w_parameters']:
                    if k.lower() == 'p3qaa_a_*':
                        cur_par_out += [k for k in p3_pars if (k not in cur_par_out) & ('p3qaa_a_' in k.lower())]
                    if k.lower() == 'p3qaa_bb_*':
                        cur_par_out += [k for k in p3_pars if (k not in cur_par_out) & ('p3qaa_bb_' in k.lower())]
                    if k.lower() == 'p3qaa_kd_*':
                        cur_par_out += [k for k in p3_pars if (k not in cur_par_out) & ('p3qaa_kd_' in k.lower())]

            ## reformat for output
            for p in cur_par_out:
                if p[6:] not in ret: ## this means we have a three band parameter
                    k = p[6:-4]
                    w = int(p[-3:])
                    if w == Bwave: wi = 0
                    if w == Gwave: wi = 1
                    if w == Rwave: wi = 2
                    par_data[p] = ret[k][:,:,wi]
                    par_atts[p] = {pk: par_attributes[pk] for pk in par_attributes}
                    par_atts[p]['ds_name'] = p
                else:
                    par_data[p] = ret[p[6:]]
                    par_atts[p] = {pk: par_attributes[pk] for pk in par_attributes}
                    par_atts[p]['ds_name'] = p
            ret = None
            p3qaa_computed = True
        ## end Pitarch 3 band QAA
        #############################

        #############################
        ## FAI
        if (cur_par == 'fai') | (cur_par == 'fai_rhot') | (cur_par == 'afai') | (cur_par == 'afai_rhot'):
            par_name = cur_par
            mask = False ## no water mask
            par_split = cur_par.split('_')
            if cur_par[0:3] == 'fai':
                par_attributes = {'algorithm':'Floating Algal Index, Hu et al. 2009', 'dataset':'rhos'}
                par_attributes['standard_name']='fai'
                par_attributes['long_name']='Floating Algal Index'
                par_attributes['units']="1"
                par_attributes['reference']='Hu et al. 2009'
                par_attributes['algorithm']=''
                ## wavelengths and max wavelength difference
                req_waves = [660, 865, 1610]
                fai_diff = [30, 30, 80]
            else:
                par_attributes = {'algorithm':'Alternative Floating Algal Index, Wang et Hu 2016', 'dataset':'rhos'}
                par_attributes['standard_name']='afai'
                par_attributes['long_name']='Alternative Floating Algal Index'
                par_attributes['units']="1"
                par_attributes['reference']='Wang et Hu 2016'
                par_attributes['algorithm']=''
                ## wavelengths and max wavelength difference
                req_waves = [660, 750, 865]
                fai_diff = [30, 30, 80]

            ## select bands
            required_datasets,req_waves_selected = [],[]
            ds_waves = [w for w in rhos_waves]
            ds_names = [ds for ds in rhos_ds]
            if 'rhot' in cur_par:
                par_attributes['dataset']='rhot'
                ds_waves = [w for w in rhot_waves]
                ds_names = [ds for ds in rhot_ds]

            ## number of products possible (for VIIRS use both I and M bands)
            fai_products = 1
            if 'VIIRS' in gem.gatts['sensor']: fai_products = 2

            ## find red and NIR bands to use
            fai_waves = []
            fai_datasets = []
            for pi in range(fai_products):
                fai_waves.append([])
                fai_datasets.append([])

                for ci, cw in enumerate(req_waves):
                    if 'VIIRS' in gem.gatts['sensor']:
                        band_prefix = {0:'I', 1:'M'}[pi]
                        if (par_attributes['standard_name']=='afai') & (pi == 0):
                            continue ## skip I band based one
                            if (ci == 1): band_prefix = 'M' ## override 750 wavelength to M band
                        ds_names_ = [ds for ii, ds in enumerate(ds_names) if ds[5] == band_prefix]
                        ds_waves_ = [ds_waves[ii] for ii, ds in enumerate(ds_names) if ds[5] == band_prefix]
                        wi, wv = ac.shared.closest_idx(ds_waves_, cw)
                        ds = ds_names_[wi]
                    else:
                        wi, wv = ac.shared.closest_idx(ds_waves, cw)
                        ds = ds_names[wi]
                    if wv > cw+fai_diff[ci]: continue
                    if wv < cw-fai_diff[ci]: continue
                    fai_waves[pi].append(wv)
                    fai_datasets[pi].append(ds)

            ## run through needed products
            for pi in range(fai_products):
                if len(fai_datasets[pi]) != len(req_waves): continue
                par_attributes['waves']=fai_waves[pi]

                ## output parameter name
                par_name = '{}'.format(cur_par)
                if 'VIIRS' in gem.gatts['sensor']:
                    par_name += '_{}'.format(fai_datasets[pi][0][5].upper())

                ## get data
                for di, cur_ds in enumerate(fai_datasets[pi]):
                    if di == 0: tmp_data = []
                    cur_data = 1.0 * gem.data(cur_ds)
                    tmp_data.append(cur_data)

                ## compute fai
                fai_sc = (float(par_attributes['waves'][1])-float(par_attributes['waves'][0]))/\
                         (float(par_attributes['waves'][2])-float(par_attributes['waves'][0]))
                #nir_prime = tmp_data[0] + (tmp_data[2]-tmp_data[0]) * fai_sc
                #par_data[par_name] = tmp_data[1] - nir_prime
                par_data[par_name] = tmp_data[1] - (1 - fai_sc) * tmp_data[0] - fai_sc * tmp_data[2]
                par_atts[par_name] = par_attributes
                tmp_data = None
        ## end FAI
        #############################

        #############################
        ## FAIT
        if (cur_par == 'fait'):
            par_name = cur_par
            mask = False ## no water mask
            par_split = par_name.split('_')
            par_attributes = {'algorithm':'Floating Algal Index Turbid Waters, Dogliotti et al. 2018', 'dataset':'rhos'}
            par_attributes['standard_name']='fait'
            par_attributes['long_name']='Floating Algal Index for Turbid Waters'
            par_attributes['units']="1"
            par_attributes['reference']='Dogliotti et al. 2018'
            par_attributes['algorithm']=''

            ## read config
            fait_cfg = ac.shared.import_config('{}/Shared/algorithms/Dogliotti/dogliotti_fait.txt'.format(ac.config['data_dir']), parse=True)
            fait_fai_threshold = float(fait_cfg['fait_fai_threshold'])
            fait_red_threshold = float(fait_cfg['fait_red_threshold'])
            fait_rgb_limit = float(fait_cfg['fait_rgb_limit'])
            fait_L_limit = float(fait_cfg['fait_L_limit'])

            if gem.gatts['sensor'] in ['L8_OLI', 'L9_OLI']:
                fait_a_threshold = float(fait_cfg['fait_a_threshold_OLI'])
            elif gem.gatts['sensor'] in ['S2A_MSI', 'S2B_MSI', 'S2C_MSI']:
                fait_a_threshold = float(fait_cfg['fait_a_threshold_MSI'])
            elif 'VIIRS' in gem.gatts['sensor']:
                fait_a_threshold = float(fait_cfg['fait_a_threshold_VIIRS'])
                fait_red_threshold = float(fait_cfg['fait_red_threshold_VIIRS'])
            else:
                print('Parameter {} not configured for {}.'.format(par_name,gem.gatts['sensor']))
                continue

            ## add to parameter attributes
            par_attributes['fai_threshold'] = fait_fai_threshold
            par_attributes['red_threshold'] = fait_red_threshold
            par_attributes['rgb_limit'] = fait_rgb_limit
            par_attributes['L_limit'] = fait_L_limit
            par_attributes['a_threshold'] = fait_a_threshold

            ## wavelengths and max wavelength difference
            fai_diff = fait_cfg['fait_diff']
            req_waves = fait_cfg['fait_wave']

            ## select bands
            required_datasets,req_waves_selected = [],[]
            ds_waves = [w for w in rhos_waves]
            ds_names = [ds for ds in rhos_ds]
            if cur_par=='fai_rhot':
                par_attributes['dataset']='rhot'
                ds_waves = [w for w in rhot_waves]
                ds_names = [ds for ds in rhot_ds]

            ## number of products possible (for VIIRS use both I and M bands)
            fai_products = 1
            if 'VIIRS' in gem.gatts['sensor']: fai_products = 2

            ## find bands to use
            fai_waves = []
            fai_datasets = []
            for pi in range(fai_products):
                fai_waves.append([])
                fai_datasets.append([])

                for ci, cw in enumerate(req_waves):
                    if 'VIIRS' in gem.gatts['sensor']:
                        band_prefix = {0:'I', 1:'M'}[pi]
                        ## override for M bands to make RGB
                        if ci <= 1: band_prefix = 'M'
                        ds_names_ = [ds for ii, ds in enumerate(ds_names) if ds[5] == band_prefix]
                        ds_waves_ = [ds_waves[ii] for ii, ds in enumerate(ds_names) if ds[5] == band_prefix]
                        wi, wv = ac.shared.closest_idx(ds_waves_, cw)
                        ds = ds_names_[wi]
                    else:
                        wi, wv = ac.shared.closest_idx(ds_waves, cw)
                        ds = ds_names[wi]
                    if wv > cw+fai_diff[ci]: continue
                    if wv < cw-fai_diff[ci]: continue
                    fai_waves[pi].append(wv)
                    fai_datasets[pi].append(ds)

            ## run through needed products
            for pi in range(fai_products):
                if len(fai_datasets[pi]) != len(req_waves): continue
                par_attributes['waves']=fai_waves[pi]

                ## output parameter name
                par_name = '{}'.format(cur_par)
                if 'VIIRS' in gem.gatts['sensor']:
                    par_name += '_{}'.format(fai_datasets[pi][-1][5].upper())

                ## get data
                for di, cur_ds in enumerate(fai_datasets[pi]):
                    if di == 0: tmp_data = []
                    cur_data = 1.0 * gem.data(cur_ds)
                    tmp_data.append(cur_data)

                ## compute fait
                fai_sc = (float(par_attributes['waves'][3])-float(par_attributes['waves'][2]))/\
                         (float(par_attributes['waves'][4])-float(par_attributes['waves'][2]))
                nir_prime = tmp_data[2] + \
                            (tmp_data[4]-tmp_data[2]) * fai_sc
                par_data[par_name] = tmp_data[3] - nir_prime
                par_atts[par_name] = par_attributes
                nir_prime = None

                ## make lab coordinates
                for i in range(3):
                    data = ac.shared.datascl(tmp_data[i], dmin=0, dmax=fait_rgb_limit)
                    if i == 0:
                        rgb = data
                    else:
                        rgb = np.dstack((rgb,data))
                    data = None
                lab = skimage.color.rgb2lab(rgb)
                rgb = None

                ## check FAI > 0
                par_data[par_name][par_data[par_name] >= fait_fai_threshold] = 1.0
                par_data[par_name][par_data[par_name] < fait_fai_threshold] = 0.0

                ## check turbidity based on red threshold
                par_data[par_name][tmp_data[2] > fait_red_threshold] = 0.0

                ## check L and a
                par_data[par_name][lab[:,:,0] >= fait_L_limit] = 0.0
                par_data[par_name][lab[:,:,1] >= fait_a_threshold] = 0.0

                ## optionally output Lab data
                if fait_cfg['fait_output_lab']:
                    for ii, lb in enumerate(['L','a', 'b']):
                        par_data['fait_Lab_{}'.format(lb)] = lab[:,:,ii]
                        par_atts['fait_Lab_{}'.format(lb)] = {}
            lab = None
        ## end FAIT
        #############################

        #############################
        ## NDVI
        if (cur_par == 'ndvi') | (cur_par.split('_')[0] == 'ndvi'):
            mask = False ## no water mask
            par_split = cur_par.split('_')
            par_attributes = {'algorithm':'NDVI', 'dataset':'rhos'}
            par_attributes['reference']=''
            par_attributes['algorithm']=''

            ## wavelengths and max wavelength difference
            ndvi_diff = [40, 100]
            req_waves = [660,865]

            ## select bands
            required_datasets,req_waves_selected = [],[]
            ds_waves = [w for w in rhos_waves]
            ds_names = [ds for ds in rhos_ds]
            if len(par_split) == 2:
                par_attributes['dataset']=par_split[1]
                ds_names = [ds for ds in gem.datasets if '{}_'.format(par_split[1]) in ds]
                ds_waves = [int(ds.split('_')[-1]) for ds in ds_names]

            ## number of products possible (for VIIRS use both I and M bands)
            ndvi_products = 1
            if 'VIIRS' in gem.gatts['sensor']: ndvi_products = 2

            ## find red and NIR bands to use
            ndvi_waves = []
            ndvi_datasets = []
            for pi in range(ndvi_products):
                ndvi_waves.append([])
                ndvi_datasets.append([])

                for ci, cw in enumerate(req_waves):
                    if 'VIIRS' in gem.gatts['sensor']:
                        band_prefix = {0:'I', 1:'M'}[pi]
                        ds_names_ = [ds for ii, ds in enumerate(ds_names) if ds[5] == band_prefix]
                        ds_waves_ = [ds_waves[ii] for ii, ds in enumerate(ds_names) if ds[5] == band_prefix]
                        wi, wv = ac.shared.closest_idx(ds_waves_, cw)
                        ds = ds_names_[wi]
                    else:
                        wi, wv = ac.shared.closest_idx(ds_waves, cw)
                        ds = ds_names[wi]
                    if wv > cw+ndvi_diff[ci]: continue
                    if wv < cw-ndvi_diff[ci]: continue
                    ndvi_waves[pi].append(wv)
                    ndvi_datasets[pi].append(ds)


            ## run through needed products
            for pi in range(ndvi_products):
                if len(ndvi_datasets[pi]) != len(req_waves): continue
                par_attributes['waves']=ndvi_waves[pi]

                ## output parameter name
                par_name = '{}'.format(cur_par)
                if 'VIIRS' in gem.gatts['sensor']:
                    par_name += '_{}'.format(ndvi_datasets[pi][0][5].upper())

                ## get data
                for di, cur_ds in enumerate(ndvi_datasets[pi]):
                    if di == 0: tmp_data = []
                    cur_data = 1.0 * gem.data(cur_ds)
                    tmp_data.append(cur_data)

                ## compute ndvi
                par_data[par_name] = (tmp_data[1]-tmp_data[0])/\
                                       (tmp_data[1]+tmp_data[0])
                par_atts[par_name] = par_attributes
                tmp_data = None
        ## end NDVI
        #############################

        #############################
        ## NDCI
        if (cur_par == 'ndci'):
            par_name = cur_par
            mask = True ## apply non water mask
            par_attributes = {'algorithm':'Mishra et al. 2014, NDCI', 'dataset':'rhos'}
            par_attributes['standard_name']='ndci'
            par_attributes['long_name']='Normalised Difference Chlorophyll Index'
            par_attributes['units']="1"
            par_attributes['reference']='Mishra et al. 2014'
            par_attributes['algorithm']=''

            required_datasets,req_waves_selected = [],[]
            ds_waves = [w for w in rhos_waves]

            ### get required datasets
            if gem.gatts['sensor'] not in ['S2A_MSI', 'S2B_MSI', 'S2C_MSI', 'S3A_OLCI', 'S3B_OLCI', 'EN1_MERIS'] + ac.config['hyper_sensors']:
                print('Parameter {} not configured for {}.'.format(par_name,gem.gatts['sensor']))
                continue

            req_waves = [670,705]
            for i, reqw in enumerate(req_waves):
                widx,selwave = ac.shared.closest_idx(ds_waves, reqw)
                if abs(float(selwave)-float(reqw)) > 10: continue
                required_datasets.append(rhos_ds[widx])
                req_waves_selected.append(selwave)
            par_attributes['waves']=req_waves_selected
            if len(required_datasets) != len(req_waves): continue
            ## get data
            for di, cur_ds in enumerate(required_datasets):
                if di == 0: tmp_data = []
                cur_data = 1.0 * gem.data(cur_ds)
                tmp_data.append(cur_data)
            ## compute ndci
            par_data[par_name] = (tmp_data[1]-tmp_data[0])/\
                                   (tmp_data[1]+tmp_data[0])
            par_atts[par_name] = par_attributes
            tmp_data = None
        ## end NDCI
        #############################

        #############################
        ## SLH
        if (cur_par == 'slh'):
            par_name = cur_par
            mask = True ## apply non water mask
            par_attributes = {'algorithm':'Kudela et al. 2015, SLH', 'dataset':'rhos'}
            par_attributes['standard_name']='slh'
            par_attributes['long_name']='Scattering Line Height'
            par_attributes['units']="1"
            par_attributes['reference']='Kudela et al. 2015'
            par_attributes['algorithm']=''
            required_datasets,req_waves_selected = [],[]
            ds_waves = [w for w in rhos_waves]

            ### get required datasets
            if gem.gatts['sensor'] not in ['S2A_MSI', 'S2B_MSI', 'S2C_MSI', 'S3A_OLCI', 'S3B_OLCI'] + ac.config['hyper_sensors']:
                print('Parameter {} not configured for {}.'.format(par_name,gem.gatts['sensor']))
                continue

            req_waves = [670,705,780]
            for i, reqw in enumerate(req_waves):
                widx,selwave = ac.shared.closest_idx(ds_waves, reqw)
                if abs(float(selwave)-float(reqw)) > 10: continue
                required_datasets.append(rhos_ds[widx])
                req_waves_selected.append(selwave)
            par_attributes['waves']=req_waves_selected
            if len(required_datasets) != len(req_waves): continue

            ## get data
            for di, cur_ds in enumerate(required_datasets):
                if di == 0: tmp_data = []
                cur_data = 1.0 * gem.data(cur_ds)
                tmp_data.append(cur_data)
            ratio = (tmp_data[2]-tmp_data[0]) / \
                    (req_waves_selected[2]+req_waves_selected[0])
            par_data[par_name] = tmp_data[1] - (tmp_data[0] + (ratio)*(req_waves_selected[1]+req_waves_selected[0]))
            par_atts[par_name] = par_attributes
            tmp_data = None
        ## end SLH
        #############################

        #############################
        ## OLH
        if (cur_par == 'olh'):
            par_name = cur_par
            mask = True ## apply non water mask

            par_attributes = {'algorithm':'Castagna et al. 2020'}
            par_attributes['standard_name']='olh'
            par_attributes['long_name']='Orange Line Height'
            par_attributes['units']="1"
            par_attributes['reference']='Castagna et al. 2020'
            par_attributes['algorithm']=''

            ### get required datasets
            if gem.gatts['sensor'] not in ['L8_OLI', 'L9_OLI', 'EO1_ALI']:
                print('Parameter {} not configured for {}.'.format(cur_par,gem.gatts['sensor']))
                continue

            if gem.gatts['sensor'] == 'L8_OLI': req_waves = [561,613,655]
            if gem.gatts['sensor'] == 'L9_OLI': req_waves = [561,613,654]
            if gem.gatts['sensor'] == 'EO1_ALI': req_waves = [561,613,655]

            #required_datasets = ['rhos_{}'.format(w) for w in req_waves]
            required_datasets = [[ds for ds in rhos_ds if ('{:.0f}'.format(w) in ds)][0] for w in req_waves]

            ## get data
            for di, cur_ds in enumerate(required_datasets):
                if di == 0: tmp_data = []
                cur_data = 1.0 * gem.data(cur_ds)
                tmp_data.append(cur_data)

            ## compute parameter
            ow = (float(req_waves[2])-req_waves[1])/(float(req_waves[2])-float(req_waves[0]))
            par_data[par_name] = tmp_data[0]*ow + tmp_data[2]*(1-ow)
            par_data[par_name] = tmp_data[1]-par_data[par_name]
            tmp_data = None
            par_atts[par_name] = par_attributes
        ## end OLH
        #############################

        #############################
        ## Hue Angle
        if (cur_par[0:9] == 'hue_angle'):
            if 'hue_angle_pitarch' == cur_par.lower():
                hue_coeff = ac.parameters.pitarch.coef_hue_angle()
                algo_ref = 'Pitarch et al., in prep.'
            else:
                hue_coeff = ac.parameters.vanderwoerd.coef_hue_angle()
                algo_ref = 'Van der Woerd et al., 2018'

            if gem.gatts['sensor'] not in hue_coeff:
                print('Parameter {} not configured for {}.'.format(cur_par, gem.gatts['sensor']))
                continue

            par_name = cur_par
            mask = True ## apply non water mask

            par_attributes = {'algorithm':'Hue Angle', 'dataset':'rhos'}
            par_attributes['standard_name']='hue_angle'
            par_attributes['long_name']='Hue Angle'
            par_attributes['units']='degrees'
            par_attributes['reference']=algo_ref
            par_attributes['algorithm']=''

            req_waves = hue_coeff[gem.gatts['sensor']]['req_waves']
            hac = hue_coeff[gem.gatts['sensor']]

            required_datasets,req_waves_selected = [],[]
            ds_waves = [w for w in rhos_waves]
            for i, reqw in enumerate(req_waves):
                widx,selwave = ac.shared.closest_idx(ds_waves, reqw)
                if abs(float(selwave)-float(reqw)) > 10: continue
                required_datasets.append(rhos_ds[widx])
                req_waves_selected.append(selwave)
            par_attributes['waves']=req_waves_selected
            if len(required_datasets) != len(req_waves): continue

            ## get data
            for di, cur_ds in enumerate(required_datasets):
                if di == 0: tmp_data = []
                cur_data = 1.0 * gem.data(cur_ds)
                tmp_data.append(cur_data)

            ## compute hue angle
            yw = 1/3.
            xw = 1/3.
            for iw, w in enumerate(req_waves_selected):
                idx, w_ = ac.shared.closest_idx(hac['lambda'], req_waves_selected[iw])
                if iw == 0:
                    X = tmp_data[iw] * hac['X'][idx]
                    Y = tmp_data[iw] * hac['Y'][idx]
                    Z = tmp_data[iw] * hac['Z'][idx]
                else:
                    X += tmp_data[iw] * hac['X'][idx]
                    Y += tmp_data[iw] * hac['Y'][idx]
                    Z += tmp_data[iw] * hac['Z'][idx]
                X[np.where(mask)] = np.nan
                Y[np.where(mask)] = np.nan
                Z[np.where(mask)] = np.nan
                den = (X+Y+Z)
                x = X/den
                y = Y/den
                den = None

                ## calculate alpha
                alpha = np.mod(np.arctan2(y-yw, x-xw),2*np.pi)
                x,y = None, None
                alpha*=(180/np.pi)
                hues_100 = alpha/100.
                corr = (hac['coef'][0] * np.power(hues_100,5)) + \
                       (hac['coef'][1] * np.power(hues_100,4)) + \
                       (hac['coef'][2] * np.power(hues_100,3)) + \
                       (hac['coef'][3] * np.power(hues_100,2)) + \
                       (hac['coef'][4] * hues_100) + hac['coef'][5]
                hues_100 = None
                alpha += corr
                corr = None
                par_data[par_name] = alpha
                par_atts[par_name] = par_attributes
                alpha = None
        ## end Hue Angle
        #############################

        ## continue if parameter not computed
        if len(par_data) == 0:
            print('Parameter {} not computed'.format(cur_par))
            continue

        ## write parameters in par_data
        for cur_ds in par_data:
            ## add mask
            if (mask) & (setu['l2w_mask_water_parameters']): par_data[cur_ds][(l2_flags & flag_value)!=0] = np.nan
            ## write to NetCDF
            if verbosity > 1: print('Writing {}'.format(cur_ds))
            gemo.write(cur_ds, par_data[cur_ds], ds_att = par_atts[cur_ds])

        ## clear data
        par_data = None
        par_atts = None
    ## end parameter loop

    ## close input NetCDF
    if close:
        gem.close()
        gem = None

    ## close output NetCDF
    gemo.close()
    gemo = None

    if verbosity>0: print('Wrote {}'.format(ofile))

    ## return file path
    return(ofile)
