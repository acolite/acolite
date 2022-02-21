## def acolite_l2w
## new L2W parameter computation for L2R generic extracted miniscene
## written by Quinten Vanhellemont, RBINS
## 2021-03-09
## modifications: 2021-12-08 (QV) added nc_projection

def acolite_l2w(gem,
                settings = None,
                sub = None,
                target_file = None,
                output = None,
                load_data = True,
                return_gem = False,
                copy_datasets = ['lon', 'lat'],
                new = True,
                verbosity=5):

    import os
    import numpy as np
    import acolite as ac
    import scipy.ndimage
    import skimage.color

    ## read gem file if NetCDF
    if type(gem) is str:
        gemf = '{}'.format(gem)
        gem = ac.gem.read(gem, sub=sub, load_data=load_data)
        if 'nc_projection' in gem:
            nc_projection = gem['nc_projection']
        else:
            nc_projection = None
    gemf = gem['gatts']['gemfile']

    ## set up output file
    if target_file is None:
        output_name = os.path.basename(gemf).replace('.nc', '')
        output_name = output_name.replace('_L2R', '_L2W')
        if 'output' in settings: output = settings['output']
        odir = output if output is not None else os.path.dirname(gemf)
        ofile = '{}/{}.nc'.format(odir, output_name)
    else:
        ofile = '{}'.format(target_file)
    gem['gatts']['ofile'] = ofile

    ## combine default and user defined settings
    setu = ac.acolite.settings.parse(gem['gatts']['sensor'], settings=settings)

    ## get rhot and rhos wavelengths
    rhot_ds = [ds for ds in gem['datasets'] if 'rhot_' in ds]
    rhot_waves = [int(ds.split('_')[-1]) for ds in rhot_ds]
    if len(rhot_waves) == 0: print('{} is probably not an ACOLITE L2R file: {} rhot datasets.'.format(gemf, len(rhot_ds)))

    rhos_ds = [ds for ds in gem['datasets'] if 'rhos_' in ds]
    rhos_waves = [int(ds.split('_')[-1]) for ds in rhos_ds]
    if len(rhos_waves) == 0: print('{} is probably not an ACOLITE L2R file: {} rhos datasets.'.format(gemf, len(rhos_ds)))

    if gem['gatts']['acolite_file_type'] != 'L2R':
        print('Only L2W processing of ACOLITE L2R files supported.')
        print('{} is a "{}" file'.format(gemf, gem['gatts']['acolite_file_type']))
        return(None)

    ## read rsr
    hyper = False
    ## hyperspectral
    if gem['gatts']['sensor'] in ac.hyper_sensors:
        hyper = True
        if gem['gatts']['sensor']=='DESIS_HSI':
            ### DESIS RSR and RSR file are version-specific
            rsrd = ac.shared.rsr_dict(f"{gem['gatts']['sensor']}_{gem['gatts']['version']}")
            # Restore sensor key without version
            rsrd[gem['gatts']['sensor']] = rsrd[f"{gem['gatts']['sensor']}_{gem['gatts']['version']}"]
            del rsrd[f"{gem['gatts']['sensor']}_{gem['gatts']['version']}"]
        else:
            rsr = ac.shared.rsr_hyper(gem['gatts']['band_waves'], gem['gatts']['band_widths'])
            rsrd = ac.shared.rsr_dict(rsrd={gem['gatts']['sensor']:{'rsr':rsr}})
    else:
        rsrd = ac.shared.rsr_dict(sensor=gem['gatts']['sensor'])

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
            if grouping_key not in gem['gatts']['auto_grouping']:
                gem['gatts']['auto_grouping'] += ':{}'.format(grouping_key)
            ## run through bands
            for b in rsrd[gem['gatts']['sensor']]['wave_name']:
                w = rsrd[gem['gatts']['sensor']]['wave_nm'][b]
                wn = rsrd[gem['gatts']['sensor']]['wave_name'][b]
                if (w < nechad_range[0]) or (w > nechad_range[1]): continue
                if wn == 'nan': continue
                setu['l2w_parameters'].append(k.replace('_*', '_{}').format(wn))

    ## compute flag value to mask for water products
    flag_value = 0
    if setu['l2w_mask']:
        flag_value = 2**setu['flag_exponent_swir']
        if setu['l2w_mask_cirrus']: flag_value += 2**setu['flag_exponent_cirrus']
        if setu['l2w_mask_high_toa']: flag_value += 2**setu['flag_exponent_toa']
        if setu['l2w_mask_negative_rhow']: flag_value += 2**setu['flag_exponent_negative']

    ## compute mask
    ## non water/swir threshold
    cidx,cwave = ac.shared.closest_idx(rhot_waves, setu['l2w_mask_wave'])
    cur_par = 'rhot_{}'.format(cwave)
    if cur_par in gem['data']:
        cur_data = 1.0 * gem['data'][cur_par]
    else:
        cur_data = ac.shared.nc_data(gemf, cur_par, sub=sub).data
    if setu['l2w_mask_smooth']:
        cur_data = ac.shared.fillnan(cur_data)
        cur_data = scipy.ndimage.gaussian_filter(cur_data, setu['l2w_mask_smooth_sigma'], mode='reflect')
    cur_mask = cur_data > setu['l2w_mask_threshold']
    cur_data = None
    l2_flags = cur_mask.astype(np.int32)*(2**setu['flag_exponent_swir'])
    cur_mask = None
    ## cirrus masking
    cidx,cwave = ac.shared.closest_idx(rhot_waves, setu['l2w_mask_cirrus_wave'])
    if np.abs(cwave - setu['l2w_mask_cirrus_wave']) < 5:
        cur_par = 'rhot_{}'.format(cwave)
        if cur_par in gem['data']:
            cur_data = 1.0 * gem['data'][cur_par]
        else:
            cur_data = ac.shared.nc_data(gemf, cur_par, sub=sub).data
        if setu['l2w_mask_smooth']:
            cur_data = ac.shared.fillnan(cur_data)
            cur_data = scipy.ndimage.gaussian_filter(cur_data, setu['l2w_mask_smooth_sigma'], mode='reflect')
        cirrus_mask = cur_data > setu['l2w_mask_cirrus_threshold']
        cirrus = None
        l2_flags += cirrus_mask.astype(np.int32)*(2**setu['flag_exponent_cirrus'])
        cirrus_mask = None
    else:
        if verbosity > 2: print('No suitable band found for cirrus masking.')
    ## TOA out of limit
    toa_mask = None
    for ci, cur_par in enumerate(rhot_ds):
        if cur_par in gem['data']:
            cur_data = 1.0 * gem['data'][cur_par]
        else:
            cur_data = ac.shared.nc_data(gemf, cur_par, sub=sub).data
        if ci == 0:
            outmask = np.isnan(cur_data)
        else:
            outmask = (outmask) | (np.isnan(cur_data))
        if setu['l2w_mask_smooth']:
            cur_data = ac.shared.fillnan(cur_data)
            cur_data = scipy.ndimage.gaussian_filter(cur_data, setu['l2w_mask_smooth_sigma'], mode='reflect')
        if toa_mask is None: toa_mask = np.zeros(cur_data.shape).astype(bool)
        toa_mask = (toa_mask) | (cur_data > setu['l2w_mask_high_toa_threshold'])
    l2_flags = (l2_flags) | (toa_mask.astype(np.int32)*(2**setu['flag_exponent_toa']))
    toa_mask = None
    l2_flags = (l2_flags) | (outmask.astype(np.int32)*(2**setu['flag_exponent_outofscene']))
    outmask = None

    ## negative rhos
    neg_mask = None
    for ci, cur_par in enumerate(rhos_ds):
        if rhos_waves[ci]<setu['l2w_mask_negative_wave_range'][0]: continue
        if rhos_waves[ci]>setu['l2w_mask_negative_wave_range'][1]: continue
        if cur_par in gem['data']:
            cur_data = 1.0 * gem['data'][cur_par]
        else:
            cur_data = ac.shared.nc_data(gemf, cur_par, sub=sub).data
        #if setu['l2w_mask_smooth']: cur_data = scipy.ndimage.gaussian_filter(cur_data, setu['l2w_mask_smooth_sigma'])
        if neg_mask is None: neg_mask = np.zeros(cur_data.shape).astype(bool)
        neg_mask = (neg_mask) | (cur_data < 0)
    l2_flags = (l2_flags) | (neg_mask.astype(np.int32)*(2**setu['flag_exponent_negative']))
    neg_mask = None

    ## list datasets to copy over from L2R
    for cur_par in gem['data']:
        if cur_par in copy_datasets: continue

        ## add rhow / Rrs from rhos
        if ('rhos_' in cur_par):
            if (('rhow_*' in setu['l2w_parameters'])) |\
                (cur_par.replace('rhos_', 'rhow_') in setu['l2w_parameters']):
                copy_datasets.append(cur_par.replace('rhos_', 'rhow_'))
            if (('Rrs_*' in setu['l2w_parameters'])) |\
                (cur_par.replace('rhos_', 'Rrs_') in setu['l2w_parameters']):
                copy_datasets.append(cur_par.replace('rhos_', 'Rrs_'))

        ## add existing par or evaluate wildcards
        if (cur_par in setu['l2w_parameters']):
            copy_datasets.append(cur_par)
        elif (('rhot_*' in setu['l2w_parameters']) & ('rhot_' in cur_par)):
            copy_datasets.append(cur_par)
        elif (('rhos_*' in setu['l2w_parameters']) & ('rhos_' in cur_par)):
            copy_datasets.append(cur_par)
        elif (('rhorc_*' in setu['l2w_parameters']) & ('rhorc_' in cur_par)):
            copy_datasets.append(cur_par)
        elif (('bt*' in setu['l2w_parameters']) & ('bt' == cur_par.lower()[0:2])):
            copy_datasets.append(cur_par)

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
        if 'rhow_' in cur_par:
            factor = 1.0
            cur_tag = cur_par.replace('rhow_','rhos_')
            mask = True
            att_add = {'algorithm':'Water reflectance', 'dataset':'rhos'}
            att_add['reference']=''
            att_add['algorithm']=''

        ## if data already read copy here
        if cur_tag in gem['data']:
            cur_data = factor * gem['data'][cur_tag]
            cur_att = gem['atts'][cur_tag]
        else:
            if cur_tag not in ac.shared.nc_datasets(gemf): continue
            cur_d, cur_att = ac.shared.nc_data(gemf, cur_tag, sub=sub, attributes=True)
            cur_data = factor * cur_d.data
            cur_data[cur_d.mask] = np.nan
            cur_d = None
        ## apply mask to Rrs and rhow
        if (mask) & (setu['l2w_mask_water_parameters']): cur_data[(l2_flags & flag_value)!=0] = np.nan
        if verbosity > 1: print('Writing {}'.format(cur_par))
        ## add attributes
        for k in att_add: cur_att[k] = att_add[k]
        ac.output.nc_write(ofile, cur_par, cur_data, dataset_attributes=cur_att,
                           attributes=gem['gatts'], new=new, nc_projection=nc_projection,
                           netcdf_compression=setu['netcdf_compression'],
                           netcdf_compression_level=setu['netcdf_compression_level'],
                           netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
        cur_data = None
        new = False

    ## write l2 flags
    ac.output.nc_write(ofile, 'l2_flags', l2_flags, attributes=gem['gatts'], new=new,
                        nc_projection=nc_projection,
                        netcdf_compression=setu['netcdf_compression'],
                        netcdf_compression_level=setu['netcdf_compression_level'])
    if return_gem: gem['data']['l2_flags'] = l2_flags
    new = False

    qaa_computed, p3qaa_computed = False, False
    ## parameter loop
    ## compute other parameters
    for cur_par in setu['l2w_parameters']:
        if cur_par.lower() in ['rhot_*', 'rhos_*', 'rrs_*', 'rhow_*', 'rhorc_*', '', ' ']: continue ## we have copied these above
        if cur_par.lower() in [ds.lower() for ds in ac.shared.nc_datasets(ofile)]: continue ## parameter already in output dataset (would not work if we are appending subsets to the ncdf)
        if cur_par.lower()[0:2] == 'bt': continue

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
            for b in rsrd[gem['gatts']['sensor']]['rsr_bands']:
                if (rsrd[gem['gatts']['sensor']]['wave_name'][b] == str(cw)): nechad_band = b

            A_Nechad, C_Nechad = None, None

            ## band center
            if nechad_parameter == 'centre':
                par_name = '{}_Nechad{}_{}'.format(nechad_par, cal_year, cw)
                nechad_dict = ac.parameters.nechad.coef_hyper(nechad_par)
                didx,algwave = ac.shared.closest_idx(nechad_dict['wave'], cw)
                A_Nechad = nechad_dict['A'][didx]
                C_Nechad = nechad_dict['C'][didx]

            ## resampled to band
            if nechad_parameter == 'resampled':
                par_name = '{}_Nechad{}Ave_{}'.format(nechad_par, cal_year, cw)
                nechad_dict = ac.parameters.nechad.coef_hyper(nechad_par)
                ## resample parameters to band
                cdct = ac.shared.rsr_convolute_dict(nechad_dict['wave']/1000, nechad_dict['C'], rsrd[gem['gatts']['sensor']]['rsr'])
                adct = ac.shared.rsr_convolute_dict(nechad_dict['wave']/1000, 1/nechad_dict['A'], rsrd[gem['gatts']['sensor']]['rsr'])
                adct = {k:1/adct[k] for k in adct}
                A_Nechad = adct[nechad_band]
                C_Nechad = cdct[nechad_band]

            ## resampled to band by Nechad 2016
            if nechad_parameter == '2016':
                par_name = '{}_Nechad2016_{}'.format(nechad_par, cw)
                par_attributes['algorithm']='2016 calibration'
                nechad_dict = ac.parameters.nechad.coef_2016()
                for k in nechad_dict:
                    if k['sensor'] != gem['gatts']['sensor'].split('_')[-1]: continue
                    if k['band'] != 'B{}'.format(nechad_band): continue
                    if k['par'] != nechad_par: continue
                    A_Nechad = k['A']
                    C_Nechad = k['C']

            ## if we have A and C we can continue
            if (A_Nechad is not None) & (C_Nechad is not None):
                cur_ds = 'rhos_{}'.format(cw)
                if cur_ds in gem['data']:
                    cur_data = 1.0 * gem['data'][cur_ds]
                else:
                    cur_data = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
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
            for b in rsrd[gem['gatts']['sensor']]['rsr_bands']:
                if (rsrd[gem['gatts']['sensor']]['wave_name'][b] == str(rw)): red_band = b
                if (rsrd[gem['gatts']['sensor']]['wave_name'][b] == str(nw)): nir_band = b

            ## store settings in atts
            for k in cfg: par_attributes[k] = cfg[k]
            par_attributes['wave_red'] = rw
            par_attributes['wave_nir'] = nw

            ## read red data
            cur_ds = 'rhos_{}'.format(rw)
            if cur_ds in gem['data']:
                red = 1.0 * gem['data'][cur_ds]
            else:
                red = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
            tur = (par_attributes['A_T_red'] * red) / (1.-red/par_attributes['C_T_red'])

            ## read nir data
            cur_ds = 'rhos_{}'.format(nw)
            if cur_ds in gem['data']:
                nir = 1.0 * gem['data'][cur_ds]
            else:
                nir = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
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
            for b in rsrd[gem['gatts']['sensor']]['rsr_bands']:
                if (rsrd[gem['gatts']['sensor']]['wave_name'][b] == str(cw)): nechad_band = b

            A_Nechad, C_Nechad = None, None
            nechad_dict = ac.parameters.dogliotti.coef_hyper()

            ## band center
            if nechad_parameter == 'center':
                par_name = '{}_Dogliotti2022_{}'.format(nechad_par, cw)
                didx,algwave = ac.shared.closest_idx(nechad_dict['wave'], cw)
                A_Nechad = nechad_dict['A_mean'][didx]
                C_Nechad = nechad_dict['C'][didx]
            elif nechad_parameter == 'resampled':
                par_name = '{}_Dogliotti2022_{}'.format(nechad_par, cw)
                # resample parameters to band
                cdct = ac.shared.rsr_convolute_dict(nechad_dict['wave']/1000, nechad_dict['C'], rsrd[gem['gatts']['sensor']]['rsr'])
                adct = ac.shared.rsr_convolute_dict(nechad_dict['wave']/1000, 1/nechad_dict['A_mean'], rsrd[gem['gatts']['sensor']]['rsr'])
                adct = {k:1/adct[k] for k in adct}
                A_Nechad = adct[nechad_band]
                C_Nechad = cdct[nechad_band]

            ## if we have A and C we can continue
            if (A_Nechad is not None) & (C_Nechad is not None):
                cur_ds = 'rhos_{}'.format(cw)
                if cur_ds in gem['data']:
                    cur_data = 1.0 * gem['data'][cur_ds]
                else:
                    cur_data = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
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
        ## CHL_OC
        if 'chl_oc' in cur_par:
            mask = True ## water parameter so apply mask
            ## load config
            chl_oc_wl_diff = 20
            cfg = ac.parameters.chl_oc.coef()
            if gem['gatts']['sensor'] not in cfg:
                print('{} not configured for {}'.format(cur_par, gem['gatts']['sensor']))
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

            if par_name not in cfg[gem['gatts']['sensor']]:
                print('{} not configured for {}'.format(par_name, gem['gatts']['sensor']))
                continue
            par_attributes['ds_name']=par_name
            chl_dct = cfg[gem['gatts']['sensor']][par_name]

            ## get bands
            blue, green = None, None
            for w in chl_dct['blue'] + chl_dct['green']:
                ci, cw = ac.shared.closest_idx(rhos_waves, int(w))
                if np.abs(cw-w) > chl_oc_wl_diff: continue
                cur_ds = 'rhos_{}'.format(cw)
                if cur_ds in gem['data']:
                    cur_data = 1.0 * gem['data'][cur_ds]
                else:
                    cur_data = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
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
        ## CHL_RE
        if (cur_par[0:6] == 'chl_re'):
            par_name = cur_par
            mask = True ## apply non water mask
            if gem['gatts']['sensor'] not in ['S2A_MSI', 'S2B_MSI', 'S3A_OLCI', 'S3B_OLCI', 'EN1_MERIS'] + ac.hyper_sensors:
                print('Parameter {} not configured for {}.'.format(par_name,gem['gatts']['sensor']))
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
                        print('Parameter {} not configured for {}.'.format(par_name,gem['gatts']['sensor']))
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
                        print('Parameter {} not configured for {}.'.format(par_name,gem['gatts']['sensor']))
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
                        selds='{}_{}'.format(par_attributes['dataset'],selwave)
                        required_datasets.append(selds)
                        req_waves_selected.append(selwave)
                    par_attributes['waves']=req_waves_selected
                    if len(required_datasets) != len(req_waves):
                        print('Required datasets not found for {}.'.format(par_name))
                        continue
                    ## get data
                    for di, cur_ds in enumerate(required_datasets):
                        if di == 0: tmp_data = []
                        if cur_ds in gem['data']:
                            cur_data  = 1.0 * gem['data'][cur_ds]
                        else:
                            cur_data  = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
                        tmp_data.append(cur_data)
                else:
                    print('Parameter {} not configured for {}.'.format(par_name,gem['gatts']['sensor']))
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
            sensor = gem['gatts']['sensor']
            if sensor not in ['L8_OLI', 'L9_OLI', 'S2A_MSI', 'S2B_MSI']:
                print('QAA not configured for {}'.format(gem['gatts']['sensor']))
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
                cur_ds = 'rhos_{}'.format(cw)
                sen_wave.append(cw)
                if ki == 0: qaa_in = {}
                if cur_ds in gem['data']:
                    cur_data  = 1.0 * gem['data'][cur_ds]
                else:
                    cur_data  = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
                ## mask data
                if (mask) & (setu['l2w_mask_water_parameters']): cur_data[(l2_flags & flag_value)!=0] = np.nan
                ## convert to Rrs
                qaa_in[k] = cur_data/np.pi

            ## get sun zenith angle
            if 'sza' in gem['data']:
                sza = gem['data']['sza']
            else:
                sza = gem['gatts']['sza']

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
            if gem['gatts']['sensor'] not in cfg:
                print('P3QAA not configured for {}'.format(gem['gatts']['sensor']))
                continue

            par_attributes = {'algorithm':'Pitarch and Vanhellemont, submitted'}

            ## read Blue Green Red data, convert to Rrs
            for k in ['B', 'G', 'R']:
                ci, cw = ac.shared.closest_idx(rhos_waves, cfg[gem['gatts']['sensor']]['center_wl'][k])
                cur_ds = 'rhos_{}'.format(cw)
                if cur_ds in gem['data']:
                    cur_data = 1.0 * gem['data'][cur_ds]
                else:
                    cur_data = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
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
            ret = ac.parameters.pitarch.p3qaa_compute(gem['gatts']['sensor'], B, G, R)

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
        if (cur_par == 'fai') | (cur_par == 'fai_rhot'):
            par_name = cur_par
            mask = False ## no water mask
            par_split = cur_par.split('_')
            par_attributes = {'algorithm':'Floating Algal Index, Hu et al. 2009', 'dataset':'rhos'}
            par_attributes['standard_name']='fai'
            par_attributes['long_name']='Floating Algal Index'
            par_attributes['units']="1"
            par_attributes['reference']='Hu et al. 2009'
            par_attributes['algorithm']=''

            ## select bands
            required_datasets,req_waves_selected = [],[]
            ds_waves = [w for w in rhos_waves]
            if cur_par=='fai_rhot':
                par_attributes['dataset']='rhot'
                ds_waves = [w for w in rhot_waves]
            ## wavelengths and max wavelength difference
            fai_diff = [10, 30, 80]
            req_waves = [660,865,1610]
            for i, reqw in enumerate(req_waves):
                widx,selwave = ac.shared.closest_idx(ds_waves, reqw)
                if abs(float(selwave)-float(reqw)) > fai_diff[i]: continue
                selds='{}_{}'.format(par_attributes['dataset'],selwave)
                required_datasets.append(selds)
                req_waves_selected.append(selwave)
            par_attributes['waves']=req_waves_selected
            if len(required_datasets) != len(req_waves): continue

            ## get data
            for di, cur_ds in enumerate(required_datasets):
                if di == 0: tmp_data = []
                if cur_ds in gem['data']:
                    cur_data  = 1.0 * gem['data'][cur_ds]
                else:
                    cur_data  = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
                tmp_data.append(cur_data)
            ## compute fai
            fai_sc = (float(par_attributes['waves'][1])-float(par_attributes['waves'][0]))/\
                     (float(par_attributes['waves'][2])-float(par_attributes['waves'][0]))
            nir_prime = tmp_data[0] + (tmp_data[2]-tmp_data[0]) * fai_sc
            par_data[par_name] = tmp_data[1] - nir_prime
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
            fait_cfg = ac.shared.import_config('{}/Shared/algorithms/Dogliotti/dogliotti_fait.txt'.format(ac.config['data_dir']))
            fait_fai_threshold = float(fait_cfg['fait_fai_threshold'])
            fait_red_threshold = float(fait_cfg['fait_red_threshold'])
            fait_rgb_limit = float(fait_cfg['fait_rgb_limit'])
            fait_L_limit = float(fait_cfg['fait_L_limit'])

            if gem['gatts']['sensor'] in ['L8_OLI', 'L9_OLI']:
                fait_a_threshold = float(fait_cfg['fait_a_threshold_OLI'])
            elif gem['gatts']['sensor'] in ['S2A_MSI', 'S2B_MSI']:
                fait_a_threshold = float(fait_cfg['fait_a_threshold_MSI'])
            else:
                print('Parameter {} not configured for {}.'.format(par_name,gem['gatts']['sensor']))
                continue

            ## add to parameter attributes
            par_attributes['fai_threshold'] = fait_fai_threshold
            par_attributes['red_threshold'] = fait_red_threshold
            par_attributes['rgb_limit'] = fait_rgb_limit
            par_attributes['L_limit'] = fait_L_limit
            par_attributes['a_threshold'] = fait_a_threshold

            ## select bands
            required_datasets,req_waves_selected = [],[]
            ds_waves = [w for w in rhos_waves]

            ## wavelengths and max wavelength difference
            fai_diff = [10, 10, 10, 30, 80]
            req_waves = [490, 560, 660, 865, 1610]
            for i, reqw in enumerate(req_waves):
                widx,selwave = ac.shared.closest_idx(ds_waves, reqw)
                if abs(float(selwave)-float(reqw)) > fai_diff[i]: continue
                selds='{}_{}'.format(par_attributes['dataset'],selwave)
                required_datasets.append(selds)
                req_waves_selected.append(selwave)
            par_attributes['waves']=req_waves_selected
            if len(required_datasets) != len(req_waves): continue

            ## get data
            for di, cur_ds in enumerate(required_datasets):
                if di == 0: tmp_data = []
                if cur_ds in gem['data']:
                    cur_data  = 1.0 * gem['data'][cur_ds]
                else:
                    cur_data  = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
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
            lab = None
        ## end FAIT
        #############################

        #############################
        ## NDVI
        if (cur_par == 'ndvi') | (cur_par == 'ndvi_rhot'):
            par_name = cur_par
            mask = False ## no water mask
            par_split = cur_par.split('_')
            par_attributes = {'algorithm':'NDVI', 'dataset':'rhos'}
            par_attributes['reference']=''
            par_attributes['algorithm']=''

            ## select bands
            required_datasets,req_waves_selected = [],[]
            ds_waves = [w for w in rhos_waves]
            if cur_par=='ndvi_rhot':
                par_attributes['dataset']='rhot'
                ds_waves = [w for w in rhot_waves]

            ## wavelengths and max wavelength difference
            ndvi_diff = [40, 65]
            req_waves = [660,865]
            for i, reqw in enumerate(req_waves):
                widx,selwave = ac.shared.closest_idx(ds_waves, reqw)
                if abs(float(selwave)-float(reqw)) > ndvi_diff[i]: continue
                selds='{}_{}'.format(par_attributes['dataset'],selwave)
                required_datasets.append(selds)
                req_waves_selected.append(selwave)
            par_attributes['waves']=req_waves_selected
            if len(required_datasets) != len(req_waves): continue
            ## get data
            for di, cur_ds in enumerate(required_datasets):
                if di == 0: tmp_data = []
                if cur_ds in gem['data']:
                    cur_data  = 1.0 * gem['data'][cur_ds]
                else:
                    cur_data  = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
                tmp_data.append(cur_data)

            ## compute ndvi
            par_data[par_name] = (tmp_data[1]-tmp_data[0])/\
                                   (tmp_data[1]+tmp_data[0])
            par_atts[par_name] = par_attributes
            tmp_data = None
        ## end NDVI

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
            if gem['gatts']['sensor'] not in ['S2A_MSI', 'S2B_MSI', 'S3A_OLCI', 'S3B_OLCI'] + ac.hyper_sensors:
                print('Parameter {} not configured for {}.'.format(par_name,gem['gatts']['sensor']))
                continue

            req_waves = [670,705]
            for i, reqw in enumerate(req_waves):
                widx,selwave = ac.shared.closest_idx(ds_waves, reqw)
                if abs(float(selwave)-float(reqw)) > 10: continue
                selds='{}_{}'.format(par_attributes['dataset'],selwave)
                required_datasets.append(selds)
                req_waves_selected.append(selwave)
            par_attributes['waves']=req_waves_selected
            if len(required_datasets) != len(req_waves): continue
            ## get data
            for di, cur_ds in enumerate(required_datasets):
                if di == 0: tmp_data = []
                if cur_ds in gem['data']:
                    cur_data  = 1.0 * gem['data'][cur_ds]
                else:
                    cur_data  = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
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
            if gem['gatts']['sensor'] not in ['S2A_MSI', 'S2B_MSI', 'S3A_OLCI', 'S3B_OLCI'] + ac.hyper_sensors:
                print('Parameter {} not configured for {}.'.format(par_name,gem['gatts']['sensor']))
                continue

            req_waves = [670,705,780]
            for i, reqw in enumerate(req_waves):
                widx,selwave = ac.shared.closest_idx(ds_waves, reqw)
                if abs(float(selwave)-float(reqw)) > 10: continue
                selds='{}_{}'.format(par_attributes['dataset'],selwave)
                required_datasets.append(selds)
                req_waves_selected.append(selwave)
            par_attributes['waves']=req_waves_selected
            if len(required_datasets) != len(req_waves): continue

            ## get data
            for di, cur_ds in enumerate(required_datasets):
                if di == 0: tmp_data = []
                if cur_ds in gem['data']:
                    cur_data  = 1.0 * gem['data'][cur_ds]
                else:
                    cur_data  = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
                tmp_data.append(cur_data)
            slh_waves = [float(ds.split('_')[1]) for ds in required_datasets]
            ratio = (tmp_data[2]-tmp_data[0]) / \
                    (slh_waves[2]+slh_waves[0])
            par_data[par_name] = tmp_data[1] - (tmp_data[0] + (ratio)*(slh_waves[1]+slh_waves[0]))
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
            if gem['gatts']['sensor'] not in ['L8_OLI', 'L9_OLI', 'EO1_ALI']:
                print('Parameter {} not configured for {}.'.format(cur_par,gem['gatts']['sensor']))
                continue

            if gem['gatts']['sensor'] == 'L8_OLI': req_waves = [561,613,655]
            if gem['gatts']['sensor'] == 'L9_OLI': req_waves = [561,613,654]
            if gem['gatts']['sensor'] == 'EO1_ALI': req_waves = [561,613,655]

            required_datasets = ['rhos_{}'.format(w) for w in req_waves]

            ## get data
            for di, cur_ds in enumerate(required_datasets):
                if di == 0: tmp_data = []
                if cur_ds in gem['data']:
                    cur_data  = 1.0 * gem['data'][cur_ds]
                else:
                    cur_data  = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
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
        if (cur_par == 'hue_angle'):
            hue_coeff = ac.parameters.vanderwoerd.coef_hue_angle()
            if gem['gatts']['sensor'] not in hue_coeff:
                print('Parameter {} not configured for {}.'.format(cur_par, gem['gatts']['sensor']))
                continue

            par_name = cur_par
            mask = True ## apply non water mask

            par_attributes = {'algorithm':'Hue Angle', 'dataset':'rhos'}
            par_attributes['standard_name']='hue_angle'
            par_attributes['long_name']='Hue Angle'
            par_attributes['units']='degrees'
            par_attributes['reference']='Van der Woerd et al., 2018'
            par_attributes['algorithm']=''


            req_waves = hue_coeff[gem['gatts']['sensor']]['req_waves']
            hac = hue_coeff[gem['gatts']['sensor']]

            required_datasets,req_waves_selected = [],[]
            ds_waves = [w for w in rhos_waves]
            for i, reqw in enumerate(req_waves):
                widx,selwave = ac.shared.closest_idx(ds_waves, reqw)
                if abs(float(selwave)-float(reqw)) > 10: continue
                selds='{}_{}'.format(par_attributes['dataset'],selwave)
                required_datasets.append(selds)
                req_waves_selected.append(selwave)
            par_attributes['waves']=req_waves_selected
            if len(required_datasets) != len(req_waves): continue

            ## get data
            for di, cur_ds in enumerate(required_datasets):
                if di == 0: tmp_data = []
                if cur_ds in gem['data']:
                    cur_data  = 1.0 * gem['data'][cur_ds]
                else:
                    cur_data  = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
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
            ac.output.nc_write(ofile, cur_ds, par_data[cur_ds], dataset_attributes=par_atts[cur_ds],
                               attributes=gem['gatts'], new=new, nc_projection=nc_projection,
                               netcdf_compression=setu['netcdf_compression'],
                               netcdf_compression_level=setu['netcdf_compression_level'],
                               netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
            ## we can also add parameter to gem
            if return_gem:
                gem['data'][cur_ds] = par_data[cur_ds]
                gem['atts'][cur_ds] = par_atts[cur_ds]

            new = False
        par_data = None
        par_atts = None
    ## end parameter loop

    ## return data or file path
    if return_gem:
        return(gem)
    else:
        return(ofile)
