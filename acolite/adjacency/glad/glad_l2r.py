## def glad_l2r
## Applies GLAD correction to ACOLITE L2R file
## GLAD developments from 2018, some results presented at ESA LPS 2019, Milan
## written by Quinten Vanhellemont, RBINS
## 2021-07-06
## modifications: 2021-10-11 (QV) added output option
##                2021-10-12 (QV) added alternative interface reflectance method
##                2021-11-09 (QV) added output of glad_x/glad_y


def glad_l2r(ncf, output = None, ofile = None,
                settings = {},
                copy_datasets = ['lat', 'lon'],
                write_rhot = True,
                glad_correction = True,
                glad_write_rhosu = True,
                glad_write_rhog = True,
                glad_write_rhoe = True,
                glad_write_xy = False,

                glad_adj_exclude_water = False,
                glad_adj_difr = False,
                glad_glint = True,
                glad_iterate = True,
                glad_start_wave = 660.,
                glad_ref1_wave = 1650.,
                glad_ref2_wave = 2200.,
                glad_neg_max = 0.01,
                glad_tau_step = 0.0025,
                glad_neg_pixel_diff = 100,

                interface_method = 'default',
                #interface_method = 'alternative',
                interface_wind = 20,
                interface_par = 'rsky_s',

                verbosity=0):

    import os
    import acolite as ac
    import numpy as np

    ## get attributes and identify sensor
    gatts = ac.shared.nc_gatts(ncf)
    nc_projection = ac.shared.nc_read_projection(ncf)
    sensor = gatts['sensor']

    sensors = ['L5_TM', 'L7_ETM', 'L8_OLI', 'L9_OLI', 'S2A_MSI','S2B_MSI']
    if sensor not in sensors:
        print('GLAD for {} not supported'.format(sensor))
        print('Supported sensors: {}'.format(', '.join(sensors)))
        return(ncf)

    ## get datasets and sensor rsr
    datasets = ac.shared.nc_datasets(ncf)
    rsrd = ac.shared.rsr_dict(sensor=sensor)

    ## get model and aot
    aot = gatts['ac_aot_550']
    lut = gatts['ac_model']

    ## get geometry and pressure
    sza, vza = gatts['sza'], gatts['vza']
    raa = gatts['raa']
    pressure = gatts['pressure']
    wind = gatts['wind']

    ## copy settings for loading LUT
    setu = {k:settings[k] for k in settings}

    ## read LUT with required parameters
    #lut_par = ['romix', 'dtott', 'utott', 'astot', 'ttot']
    #lutdw = ac.aerlut.import_luts(add_rsky=True, #, add_dutott=True, pressures = setu['luts_pressures'],
    #                              sensor=sensor, base_luts = [lut], lut_par = lut_par)
    lutdw = ac.aerlut.import_luts(add_rsky=True, par='romix+rsky_t',
                                    sensor=sensor, base_luts=[lut], pressures = setu['luts_pressures'],
                                    reduce_dimensions=setu['luts_reduce_dimensions'])

    ## get gas transmittance
    tg_dict = ac.ac.gas_transmittance(sza, vza, uoz=gatts['uoz'], uwv=gatts['uwv'], rsr=rsrd[sensor]['rsr'])

    ## make bands dataset
    bands = {}
    for bi, b in enumerate(rsrd[sensor]['rsr_bands']):
        if b not in bands:
            bands[b] = {k:rsrd[sensor][k][b] for k in ['wave_mu', 'wave_nm', 'wave_name'] if b in rsrd[sensor][k]}
            bands[b]['rhot_ds'] = 'rhot_{}'.format(bands[b]['wave_name'])
            bands[b]['rhos_ds'] = 'rhos_{}'.format(bands[b]['wave_name'])
            for k in tg_dict:
                if k not in ['wave']: bands[b][k] = tg_dict[k][b]
            bands[b]['wavelength']=bands[b]['wave_nm']
    ## end bands dataset

    ## retrieve required atmospheric parameters
    atm = {'tt_gas': tg_dict['tt_gas']}
    for par in lut_par:
        atm[par] = {b: float(lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd'][par], raa, vza, sza, wind, aot))) \
                    for b in rsrd[sensor]['rsr_bands']}

    ## set up output attributes
    gatts_out = {g:gatts[g] for g in gatts}
    gatts_out['glad_applied'] = 1

    gatts_out['glad_write_rhosu'] = 1 if glad_write_rhosu else 0
    gatts_out['glad_write_rhog'] = 1 if glad_write_rhog else 0
    gatts_out['glad_write_rhoe'] = 1 if glad_write_rhoe else 0
    gatts_out['glad_write_xy'] = 1 if glad_write_xy else 0
    gatts_out['glad_adj_exclude_water'] = 1 if glad_adj_exclude_water else 0
    gatts_out['glad_adj_difr'] = 1 if glad_adj_difr else 0
    gatts_out['glad_glint'] = 1 if glad_glint else 0
    gatts_out['glad_iterate'] = 1 if glad_iterate else 0

    gatts_out['glad_start_wave'] = glad_start_wave
    gatts_out['glad_ref1_wave'] = glad_ref1_wave
    gatts_out['glad_ref2_wave'] = glad_ref2_wave
    gatts_out['glad_neg_max'] = glad_neg_max
    gatts_out['glad_tau_step'] = glad_tau_step
    gatts_out['glad_neg_pixel_diff'] = glad_neg_pixel_diff

    if 'rhog' not in gatts_out['auto_grouping']:
        gatts_out['auto_grouping'] += ':rhog'
    if 'rhoe' not in gatts_out['auto_grouping']:
        gatts_out['auto_grouping'] += ':rhoe'
    if 'rhosu' not in gatts_out['auto_grouping']:
        gatts_out['auto_grouping'] += ':rhosu'

    ## outputfile
    if ofile is None:
        ofile = ncf.replace('_L2R.nc', '_L2R_GLAD.nc')
    if output is not None:
        ofile = '{}/{}'.format(output, os.path.basename(ofile))

    update_attributes, new = True, True
    for ds in copy_datasets:
        if ds not in datasets: continue
        if (nc_projection is not None) & (ds in ['x', 'y']): continue #gatts['projection_key']
        d, a = ac.shared.nc_data(ncf, ds, attributes=True)
        ac.output.nc_write(ofile, ds, d, dataset_attributes = a,
                           nc_projection = nc_projection,
                           attributes=gatts_out, update_attributes = update_attributes, new=new)
        update_attributes, new = False, False

    ## new implementataion of glad correction
    if glad_correction:
        dtor = np.pi / 180

        cos2omega = np.cos(sza*dtor)*np.cos(vza*dtor) + np.sin(sza*dtor)*np.sin(vza*dtor)*np.cos(raa*dtor)
        omega = np.arccos(np.sqrt(cos2omega))
        omega = np.arccos(cos2omega)/2

        ## read and resample refractive index
        refri = ac.ac.refri()
        refri_sen = ac.shared.rsr_convolute_dict(refri['wave']/1000, refri['n'], rsrd[sensor]['rsr'])

        ## compute fresnel reflectance for each n
        Rf_sen = {b: ac.ac.sky_refl(omega, n_w=refri_sen[b]) for b in refri_sen}

        ## identify bands to use
        glad_waves = [bands[b]['wave_nm'] for b in bands]
        glad_bands = [b for b in bands]
        glad_ref1_idx, glad_ref1_wv = ac.shared.closest_idx(glad_waves, glad_ref1_wave)
        glad_ref2_idx, glad_ref2_wv = ac.shared.closest_idx(glad_waves, glad_ref2_wave)
        glad_ref1_band = glad_bands[glad_ref1_idx]
        glad_ref2_band = glad_bands[glad_ref2_idx]

        ## identify start band (start in red, with usually lowest water reflectance)
        glad_start_idx, glad_start_wv = ac.shared.closest_idx(glad_waves, glad_start_wave)
        glad_order = [glad_bands[glad_start_idx]]
        glad_order += [b for b in glad_bands if b not in glad_order]

        ## load SWIR bands
        glad_ref_rhot1, rhot1_att = ac.shared.nc_data(ncf, bands[glad_ref1_band]['rhot_ds'], attributes=True)
        glad_ref_rhot2, rhot2_att = ac.shared.nc_data(ncf, bands[glad_ref2_band]['rhot_ds'], attributes=True)

        ## mask of non water
        glad_mask = glad_ref_rhot1 > 0.1
        nvalid = len(np.where(glad_mask==0)[0])
        if nvalid == 0:
            if verbosity > 0: print('No valid pixels found - skipping GLAD')
            glad_loop = False
            glad_skip = True

        ## initial conditions, use DSF aot
        caot = 1. * aot
        neg_absolute = None
        glad_iter_max = max((1,int(caot / glad_tau_step)))
        glad_skip = False

        neg_fraction = 1.0
        dneg = 101
        glad_iter = 0

        ## get alt surface reflectance
        rgi_pos = [pressure, 0, raa, vza, sza, interface_wind, np.nanmin((0.1, aot))]
        rgi_pos[1] = lutdw[lut]['ipd'][interface_par]

        int_cur = {b: float(lutdw[lut]['rgi'][b](rgi_pos)) for b in rsrd[sensor]['rsr_bands']}
        int_cur_b1 = {b: int_cur[b]/int_cur[glad_ref1_band] if int_cur[glad_ref1_band]>0 else 0 for b in rsrd[sensor]['rsr_bands'] }
        int_cur_b2 = {b: int_cur[b]/int_cur[glad_ref2_band] if int_cur[glad_ref2_band]>0 else 0 for b in rsrd[sensor]['rsr_bands']}

        ## run through bands
        for bi,b in enumerate(glad_order):
            if glad_skip: continue
            if bands[b]['rhos_ds'] not in datasets: continue
            wave = bands[b]['wave_nm']
            if verbosity > 0: print('Performing GLAD correction for band {} ({:.0f} nm)'.format(b, wave))

            ## load rhot for current band
            if b == glad_ref1_band:
                cur_rhot = 1.0 * glad_ref_rhot1
                cur_att = {a:rhot1_att[a] for a in rhot1_att}
            elif b == glad_ref2_band:
                cur_rhot = 1.0 * glad_ref_rhot2
                cur_att = {a:rhot2_att[a] for a in rhot2_att}
            else:
                cur_rhot, cur_att = ac.shared.nc_data(ncf, bands[b]['rhot_ds'], attributes=True)

            ## convert from masked array
            cur_toa_mask = cur_rhot.mask
            cur_rhot = cur_rhot.data
            cur_rhot[cur_toa_mask] = np.nan

            ## glad loop to determine best aot and relative glint/adjacency contribution
            glad_loop = True
            while glad_loop:
                ## first band
                if bi == 0:
                    ## update atm parameters if new aot is different
                    if caot != aot:
                        for par in lut_par:
                            atm[par] = {b: float(lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd'][par], raa, vza, sza, wind, caot))) \
                                        for b in rsrd[sensor]['rsr_bands']}

                        ## get alt surface reflectance
                        rgi_pos = [pressure, 0, raa, vza, sza, interface_wind, np.nanmin((0.1, aot))]
                        rgi_pos[1] = lutdw[lut]['ipd'][interface_par]

                        int_cur = {b: float(lutdw[lut]['rgi'][b](rgi_pos)) for b in rsrd[sensor]['rsr_bands']}
                        #int_cur_b1 = {b: int_cur[b]/int_cur[glad_ref1_band] for b in rsrd[sensor]['rsr_bands']}
                        #int_cur_b2 = {b: int_cur[b]/int_cur[glad_ref2_band] for b in rsrd[sensor]['rsr_bands']}
                        int_cur_b1 = {b: int_cur[b]/int_cur[glad_ref1_band] if int_cur[glad_ref1_band]>0 else 0 for b in rsrd[sensor]['rsr_bands'] }
                        int_cur_b2 = {b: int_cur[b]/int_cur[glad_ref2_band] if int_cur[glad_ref2_band]>0 else 0 for b in rsrd[sensor]['rsr_bands']}

                    ## total transmittance dtott_s utott_s
                    glad_dict = {'Tu':{}, 'Td':{}, 'T':{},
                                 'tdu':{}, 'tdd':{}, 'u_diftodir':{}, 'd_diftodir':{},
                                 'Rf_ref1':{}, 'Rf_ref2':{}, 'gc_ref1':{}, 'gc_ref2':{}}

                    ## compute for new ctau
                    for bi2,b2 in enumerate(glad_bands):
                        ## direct up and down transmittances
                        glad_dict['Tu'][b2] = np.exp(-1.*(atm['ttot'][b2]/np.cos(vza*dtor)))
                        glad_dict['Td'][b2] = np.exp(-1.*(atm['ttot'][b2]/np.cos(sza*dtor)))

                        ## two way direct transmittance
                        glad_dict['T'][b2]  = glad_dict['Tu'][b2] * glad_dict['Td'][b2]

                        ## diffuse up and down transmittances
                        glad_dict['tdu'][b2] = atm['utott'][b2]-glad_dict['Tu'][b2]
                        glad_dict['tdd'][b2] = atm['dtott'][b2]-glad_dict['Td'][b2]

                        ## diffuse to direct transmittance ratios
                        glad_dict['u_diftodir'][b2] = glad_dict['tdu'][b2] / glad_dict['Tu'][b2]
                        glad_dict['d_diftodir'][b2] = glad_dict['tdd'][b2] / glad_dict['Td'][b2]

                        ## fresnel reflectance ratio for SWIR1 and SWIR2
                        glad_dict['Rf_ref1'][b2]  = Rf_sen[b2]/Rf_sen[glad_ref1_band]
                        glad_dict['Rf_ref2'][b2]  = Rf_sen[b2]/Rf_sen[glad_ref2_band]

                        ## glint correction factor for SWIR1 and SWIR2
                        glad_dict['gc_ref1'][b2]  = glad_dict['T'][b2] * glad_dict['Rf_ref1'][b2]
                        glad_dict['gc_ref2'][b2]  = glad_dict['T'][b2] * glad_dict['Rf_ref2'][b2]

                    for bi2,b2 in enumerate(glad_bands):
                        glad_dict['gc_ref1'][b2] /= glad_dict['T'][glad_ref1_band]
                        glad_dict['gc_ref2'][b2] /= glad_dict['T'][glad_ref2_band]

                    ## do atmospheric correction
                    ## ref band 1
                    glad_ref_rhos1 = (glad_ref_rhot1/ bands[glad_ref1_band]['tt_gas']) - atm['romix'][glad_ref1_band]
                    glad_ref_rhos1 = (glad_ref_rhos1) / ((atm['dtott'][glad_ref1_band] * atm['utott'][glad_ref1_band]) + atm['astot'][glad_ref1_band]*glad_ref_rhos1)
                    glad_ref_rhos1[glad_ref_rhos1.mask] = np.nan
                    glad_ref_rhos1 = glad_ref_rhos1.data
                    ## ref band 2
                    glad_ref_rhos2 = (glad_ref_rhot2/ bands[glad_ref2_band]['tt_gas']) - atm['romix'][glad_ref2_band]
                    glad_ref_rhos2 = (glad_ref_rhos2) / ((atm['dtott'][glad_ref2_band] * atm['utott'][glad_ref2_band]) + atm['astot'][glad_ref2_band]*glad_ref_rhos2)
                    glad_ref_rhos2[glad_ref_rhos2.mask] = np.nan
                    glad_ref_rhos2 = glad_ref_rhos2.data

                    ## mask of positive SWIR reflectances
                    glad_positives = (glad_ref_rhos1>0) & (glad_ref_rhos2>0)
                    if neg_absolute is None: neg_absolute = nvalid * 1.0

                    if glad_adj_exclude_water:
                        glad_ref_rhos1_ave = np.nanmean(glad_ref_rhos1[glad_mask==1])
                        glad_ref_rhos2_ave = np.nanmean(glad_ref_rhos2[glad_mask==1])
                    else:
                        glad_ref_rhos1_ave = np.nanmean(glad_ref_rhos1)
                        glad_ref_rhos2_ave = np.nanmean(glad_ref_rhos2)

                    if glad_adj_difr:
                        ## ratio of diffuse to direct ratios
                        glad_ref_rhoe1 = glad_dict['u_diftodir'][glad_ref1_band]*(glad_ref_rhos1_ave)
                        glad_ref_rhoe2 = glad_dict['u_diftodir'][glad_ref2_band]*(glad_ref_rhos2_ave)
                    else:
                        glad_ref_rhoe1 = 1 * glad_ref_rhos1_ave
                        glad_ref_rhoe2 = 1 * glad_ref_rhos2_ave

                    #glad glint ratio
                    glad_gr = glad_dict['gc_ref1'][glad_ref2_band] # unused here

                    ## glad environment ratio
                    glad_er = glad_ref_rhoe2 / glad_ref_rhoe1

                    if interface_method == 'default':
                        glad_x = (glad_ref_rhos2 - (glad_er * glad_ref_rhos1)) /\
                                 (glad_dict['gc_ref1'][glad_ref2_band] - glad_er)
                    else:
                        #print(glad_dict['gc_ref1'])
                        #print(int_cur_b2)
                        glad_x = (glad_ref_rhos2 - (glad_er * glad_ref_rhos1)) /\
                                 (int_cur_b1[glad_ref2_band] - glad_er)

                    glad_x[glad_x<0]=0 ## do not allow negative glint
                    if not glad_glint: glad_x[:]=0 ## set to fully ignore glint contribution
                    glad_x = np.minimum(glad_x, glad_ref_rhos1)

                    ## compute the environment contribution to the SWIR
                    glad_y = (glad_ref_rhos1 - glad_x) /(glad_ref_rhoe1)
                    glad_y[glad_y<0]=0

                    ## mask negative SWIR observations
                    glad_x[glad_positives==0] = 0
                    glad_y[glad_positives==0] = 0

                if glad_skip:
                    glad_loop = False
                    continue


                ## a/c in current band
                cur_rhos = (cur_rhot/ bands[b]['tt_gas']) - atm['romix'][b]
                cur_rhos = (cur_rhos) / ((atm['dtott'][b] * atm['utott'][b]) + atm['astot'][b]*cur_rhos)
                if glad_write_rhosu: cur_rhosu = cur_rhos * 1.0

                if glad_adj_exclude_water:
                    cur_rhos_ave = np.nanmean(cur_rhos[glad_mask==1])
                else:
                    cur_rhos_ave = np.nanmean(cur_rhos)

                ## compute glint in this band
                if interface_method == 'default':
                    cur_rhog = glad_dict['gc_ref1'][b] * glad_x
                else:
                    cur_rhog = int_cur_b1[b] * glad_x

                cur_rhog[glad_mask == 1] = np.nan
                cur_rhog[cur_toa_mask] = np.nan

                ## compute environment reflectance in this band
                if glad_adj_difr:
                    cur_rhoe_w = glad_y * (glad_dict['u_diftodir'][btag])
                else:
                    cur_rhoe_w = glad_y
                cur_rhoe = cur_rhoe_w*cur_rhos_ave
                cur_rhoe[glad_mask == 1] = np.nan
                cur_rhoe[cur_toa_mask] = np.nan

                cur_rhos[glad_mask == 0] -= cur_rhog[glad_mask == 0]
                cur_rhos[glad_mask == 0] -= cur_rhoe[glad_mask == 0]

                glad_iter+=1
                if bi == 0:
                    update_attributes = False
                    nneg = len(np.where(cur_rhos[glad_mask == 0] < 0)[0])
                    if nvalid == 0:
                        glad_loop = False
                    else:
                        neg_fraction = nneg/nvalid
                        dneg = neg_absolute - nneg
                        #print(dneg, nvalid*0.01, neg_fraction)

                        glad_loop = (neg_fraction > glad_neg_max) & (glad_iter < glad_iter_max)

                        ## stop if less than 1% of pixels were affected
                        if dneg < nvalid*0.01: glad_loop = False

                        if not glad_iterate: glad_loop = False

                        if verbosity > 3: print(nneg, int(dneg), 'neg {:.1f}%'.format(neg_fraction*100), 'aot {:.3f}'.format(caot))
                        if glad_loop:
                            caot -= glad_tau_step
                            neg_absolute = 1 * nneg
                        else:
                            if glad_iterate:
                                if verbosity > 2: print('Finished band {:.0f} after {} iterations'.format(wave, glad_iter))
                                if verbosity > 2: print('Starting aot {:.3f}, end aot {:.3f}'.format(aot, caot))
                            gatts_out['ac_aot_550_dsf'] = gatts_out['ac_aot_550'] * 1.0
                            gatts_out['ac_aot_550'] = caot
                            update_attributes = True
                else:
                    glad_loop = False

                ## end glad loop
                ################
                if glad_skip: continue


            ## write current rhot
            if write_rhot:
                ac.output.nc_write(ofile, bands[b]['rhot_ds'], cur_rhot, dataset_attributes = cur_att,
                       nc_projection=nc_projection, attributes=gatts_out, update_attributes = update_attributes, new=new)
                update_attributes, new = False, False

            ## update attributes
            for k in atm: cur_att[k] = atm[k][b]
            ac.output.nc_write(ofile, bands[b]['rhos_ds'], cur_rhos, dataset_attributes = cur_att,
                                nc_projection=nc_projection, attributes=gatts_out, update_attributes = update_attributes, new=new)
            update_attributes, new = False, False


            ## write uncorrected rhos
            if glad_write_rhosu:
                ac.output.nc_write(ofile, bands[b]['rhos_ds'].replace('rhos_', 'rhosu_'), cur_rhosu,
                                  dataset_attributes = cur_att)

            ## write rhog
            if glad_write_rhog:
                ac.output.nc_write(ofile, bands[b]['rhos_ds'].replace('rhos_', 'rhog_'), cur_rhog,
                                  dataset_attributes = cur_att)

            ## write rhoe
            if glad_write_rhoe:
                ac.output.nc_write(ofile, bands[b]['rhos_ds'].replace('rhos_', 'rhoe_'), cur_rhoe,
                                  dataset_attributes = cur_att)

    ## no valid pixels found
    if glad_skip:
        if os.path.exists(ofile): os.remove(ofile)
        return(ncf)

    ## write glad_x and glad_y
    if glad_write_xy:
        ac.output.nc_write(ofile, 'glad_x', glad_x)
        ac.output.nc_write(ofile, 'glad_y', glad_y)

    ## recompute orange band
        if (sensor in ['L8_OLI', 'L9_OLI']):
            if verbosity > 2: print('Recomputing orange band')
            ## load orange band configuration
            if sensor == 'L8_OLI':
                ob_cfg = ac.shared.import_config(ac.config['data_dir']+'/L8/oli_orange.cfg')
                sensor_o = 'L8_OLI_ORANGE'
            if sensor == 'L9_OLI':
                ob_cfg = ac.shared.import_config(ac.config['data_dir']+'/L9/oli_orange.cfg')
                sensor_o = 'L9_OLI_ORANGE'

        ## load orange band data
        rsrd_o = ac.shared.rsr_dict(sensor_o)[sensor_o]
        ob = {k:rsrd_o[k]['O'] for k in ['wave_mu', 'wave_nm', 'wave_name']}
        ds = 'rhos_{}'.format(ob['wave_name'])

        d, cur_att = ac.shared.nc_data(ncf, ds, attributes=True)

        ## write uncorrected data
        if glad_write_rhosu:
            ac.output.nc_write(ofile, ds.replace('rhos_', 'rhosu_'), d, dataset_attributes=cur_att)

        ## compute orange band
        ob_data =  ac.shared.nc_data(ofile, bands['8']['rhos_ds']) * float(ob_cfg['pf'])
        ob_data += ac.shared.nc_data(ofile, bands['3']['rhos_ds']) * float(ob_cfg['gf'])
        ob_data += ac.shared.nc_data(ofile, bands['4']['rhos_ds']) * float(ob_cfg['rf'])
        ac.output.nc_write(ofile, ds, ob_data, dataset_attributes=cur_att)

    return(ofile)
