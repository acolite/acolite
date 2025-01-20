## separated glint correction function
## currently only supports one aot/model combination
## written by Quinten Vanhellemont, RBINS
## 2024-05-15
## modifications: 2024-05-16 (QV) added function, check if file exists, added new_file option
##                2024-05-22 (QV) added write keyword, update to use gem.write_ds
##                2024-07-04 (QV) added update dataset info
##                2024-07-29 (QV) fixes when gem is passed

def default(gem, settings = None, lutdw = None, write = True, new_file = False):
    import acolite as ac
    import numpy as np
    import os, time, shutil

    ## check if file path is given
    opened = False
    if type(gem) is str:
        gem = ac.gem.gem(gem)
        opened = True
    else:
        if gem.file is not None:
            gem.setup() ## update dataset info
    gemf = gem.file

    if gemf is None: write = False ## no writing if file is None

    if (write):
        if (not os.path.exists(gemf)):
            print('File does not exist {}'.format(gemf))
            if opened: gem.close()
            return

    if 'glint_correction' in gem.gatts:
        print('Glint correction already applied to {}'.format(gemf))
        if opened: gem.close()
        return

    ## create new file if requested
    ## otherwise L2R file will be updated if write = True
    if (write) & (new_file):
        gemi = '{}'.format(gemf)
        oname = os.path.basename(gemi).replace('L2R', 'L2R_GC')
        odir = os.path.dirname(gemi)
        if type(settings) is dict:
            if 'output' in settings: odir = settings['output']
        gem.close()

        ## copy file
        gemf = '{}/{}'.format(odir, oname)

        print('Creating new file: {}'.format(gemf))
        shutil.copy(gemi, gemf)
        gem = ac.gem.gem(gemf)
        opened = True

    ## combine default and user defined settings
    if settings is not None:
        ac.settings['user'] = ac.acolite.settings.parse(None, settings=settings, merge=False)
        for k in ac.settings['user']: ac.settings['run'][k] = ac.settings['user'][k]
    setu = ac.acolite.settings.parse(gem.gatts['sensor'], settings=ac.settings['user'])
    for k in setu: ac.settings['run'][k] = setu[k]  ## update run settings with user settings and sensor defaults

    ## get aerosol information, currently only fixed
    aot = gem.gatts['ac_aot_550']
    model = gem.gatts['ac_model']

    if 'aot_550' in gem.datasets:
        print('Per-pixel aot found, not yet implemented.')
        if opened:
            gem.close()
            del gem
        return

    if ('sza' in gem.datasets) & ('vza' in gem.datasets) & \
       ((('saa' in gem.datasets) & ('vaa' in gem.datasets)) | ('raa' in gem.datasets)):
        print('Per-pixel geometry found, not yet implemented.')
        print('Proceeding with scene average geometry values.')

    ## get bands dataset
    sensor, rsrd, bands = ac.gem.bands(gem)

    if 'radcor_version' in gem.gatts:
        base_datasets = ['rhos', 'rhosu']
    else:
        base_datasets = ['rhos']

    ## get parameters and load lut
    romix_par = 'romix'
    if setu['dsf_interface_reflectance']: romix_par = 'romix+rsky_t'
    add_rsky = romix_par == 'romix+rsky_t'

    ## import lut if not passed to function
    if lutdw is None:
        lutdw = ac.aerlut.import_luts(add_rsky = add_rsky, par = romix_par,
                                      sensor = sensor, lut_par = ['ttot'], return_lut_array = True)
    luts = list(lutdw.keys())

    if 'ac_lut' in gem.gatts:
        lut = gem.gatts['ac_lut']
    else:
        aer_nm = {'C': 'MOD1', 'M': 'MOD2', 'U': 'MOD3'}
        lut = [lut for lut in luts if aer_nm[model] in lut][0]

    pressure = setu['pressure_default'] # ac.settings['run']['pressure_default']
    wind = setu['wind_default'] # ac.settings['run']['wind_default']
    if 'pressure' in gem.gatts: pressure = gem.gatts['pressure']
    if 'wind' in gem.gatts: wind = gem.gatts['wind']

    for base in base_datasets:
        ## get average geometry
        raa = gem.gatts['raa']
        vza = gem.gatts['vza']
        sza = gem.gatts['sza']

        if add_rsky:
            xi = [pressure, raa, vza, sza, wind]
        else:
            xi = [pressure, raa, vza, sza]

        ## find bands for glint correction
        gc_swir1, gc_swir2 = None, None
        gc_swir1_b, gc_swir2_b = None, None
        swir1d, swir2d = 1000, 1000
        gc_user, gc_mask = None, None
        gc_user_b, gc_mask_b = None, None
        userd, maskd = 1000, 1000

        for b in bands:
            ## swir1
            sd = np.abs(bands[b]['wave_nm'] - 1600)
            if sd < 100:
                if sd < swir1d:
                    gc_swir1 = bands[b]['{}_ds'.format(base)]
                    swir1d = sd
                    gc_swir1_b = b
            ## swir2
            sd = np.abs(bands[b]['wave_nm'] - 2200)
            if sd < 100:
                if sd < swir2d:
                    gc_swir2 = bands[b]['{}_ds'.format(base)]
                    swir2d = sd
                    gc_swir2_b = b
            ## mask band
            sd = np.abs(bands[b]['wave_nm'] - setu['glint_mask_rhos_wave'])
            if sd < 100:
                if sd < maskd:
                    gc_mask = bands[b]['{}_ds'.format(base)]
                    maskd = sd
                    gc_mask_b = b
            ## user band
            if setu['glint_force_band'] is not None:
                sd = np.abs(bands[b]['wave_nm'] - setu['glint_force_band'])
                if sd < 100:
                    if sd < userd:
                        gc_user = bands[b]['{}_ds'.format(base)]
                        userd = sd
                        gc_user_b = b

            ## use user selected  band
            if gc_user is not None:
                gc_swir1, gc_swir1_b = None, None
                gc_swir2, gc_swir2_b = None, None

        ## compute total optical thickness based on aot and model
        ## currently only fixed values
        ttot_all = {}
        for b in bands:
            if len(xi) == 4:
                ttot_all[b] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd']['ttot'], xi[1], xi[2], xi[3], aot))
            if len(xi) == 5:
                ttot_all[b] = lutdw[lut]['rgi'][b]((xi[0], lutdw[lut]['ipd']['ttot'], xi[1], xi[2], xi[3], xi[4], aot))

        ## start glint correction
        if ((gc_swir1 is not None) and (gc_swir2 is not None)) or (gc_user is not None):
            t0 = time.time()
            if setu['verbosity'] > 1: print('Starting glint correction with base dataset {}'.format(base))

            ## compute scattering angle
            sza = np.radians(sza)
            vza = np.radians(vza)
            raa = np.radians(raa)

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
            gc_mask_data = gem.data(gc_mask)

            if gc_mask_data is None: ## reference rhos dataset can be missing for night time images (tgas computation goes negative)
                if setu['verbosity'] > 1: print('No glint mask could be determined.')
            else:
                sub_gc = np.where(np.isfinite(gc_mask_data) & \
                                  (gc_mask_data<=setu['glint_mask_rhos_threshold']))
                del gc_mask_data

                ## get reference bands transmittance
                for ib, b in enumerate(bands):
                    rhos_ds = bands[b]['{}_ds'.format(base)]
                    if rhos_ds not in [gc_swir1, gc_swir2, gc_user]: continue
                    if rhos_ds not in gem.datasets: continue
                    ttot_all_b = ttot_all[b] * 1.0
                    T_cur  = np.exp(-1.*(ttot_all_b/muv)) * np.exp(-1.*(ttot_all_b/mus))
                    del ttot_all_b

                    ## subset if 2d
                    T_cur_sub = T_cur[sub_gc] if len(np.atleast_2d(T_cur)) > 1 else T_cur * 1.0

                    if rhos_ds == gc_user:
                        T_USER = T_cur_sub * 1.0
                    else:
                        if rhos_ds == gc_swir1: T_SWIR1 = T_cur_sub * 1.0
                        if rhos_ds == gc_swir2: T_SWIR2 = T_cur_sub * 1.0
                    del T_cur, T_cur_sub

                ## swir band choice is made for first band
                gc_choice = False
                ## glint correction per band
                for ib, b in enumerate(bands):
                    rhos_ds = bands[b]['{}_ds'.format(base)]
                    if rhos_ds not in gem.datasets: continue
                    if b not in ttot_all: continue
                    if setu['verbosity'] > 5: print('Performing glint correction for band {} ({} nm)'.format(b, bands[b]['wave_name']))
                    ## load rhos dataset
                    cur_data, cur_att = gem.data(rhos_ds, attributes = True)

                    ## two-way direct transmittance
                    ttot_all_b = ttot_all[b] * 1.0
                    T_cur  = np.exp(-1.*(ttot_all_b/muv)) * np.exp(-1.*(ttot_all_b/mus))
                    del ttot_all_b

                    ## subset if 2d
                    T_cur_sub = T_cur[sub_gc] if len(np.atleast_2d(T_cur)) > 1 else T_cur * 1.0
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
                            swir1_rhos = gem.data(gc_swir1)[sub_gc]
                            swir2_rhos = gem.data(gc_swir2)[sub_gc]
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
                            rhog_ref = gem.data(gc_user)[sub_gc]
                            ## set negatives to 0
                            rhog_ref[rhog_ref<0] = 0

                        ## write reference glint
                        if setu['glint_write_rhog_ref']:
                            ## set x/y dims here - should have been done elsewhere
                            gem.ydim = cur_data.shape[0]
                            gem.xdim = cur_data.shape[1]
                            tmp = np.zeros((gem.ydim, gem.xdim), dtype=np.float32) + np.nan
                            tmp[sub_gc] = rhog_ref
                            ## add to gem
                            ds = 'rhog_ref'
                            if base == 'rhosu': ds = 'rhogu_ref'
                            gem.data_mem[ds] = tmp
                            gem.data_att[ds] = {}
                            ## write and clear
                            if write: gem.write_ds(ds, clear = True)
                            ## old write option
                            #gem.write('rhog_ref', tmp)
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
                    for a in bands[b]: cur_att[a] = bands[b][a]

                    ## add to gem
                    ds = rhos_ds
                    gem.data_mem[ds] = cur_data
                    gem.data_att[ds] = cur_att
                    ## write and clear
                    if write: gem.write_ds(ds, clear = True)
                    ## old write option
                    #gem.write(rhos_ds, cur_data, ds_att = cur_att)

                    ## write band glint
                    if setu['glint_write_rhog_all']:
                        ## set x/y dims here - should have been done elsewhere
                        gem.ydim = cur_data.shape[0]
                        gem.xdim = cur_data.shape[1]
                        tmp = np.zeros((gem.ydim, gem.xdim), dtype=np.float32) + np.nan
                        tmp[sub_gc] = cur_rhog
                        ## add to gem
                        ds = rhos_ds.replace('rhos_', 'rhog_')
                        if base == 'rhosu': ds = rhos_ds.replace('rhosu_', 'rhogu_')
                        gem.data_mem[ds] = tmp
                        gem.data_att[ds] = cur_att
                        ## write and clear
                        if write: gem.write_ds(ds, clear = True)
                        ## old write option
                        #gem.write(ds, tmp, ds_att=cur_att)
                        del tmp
                    del cur_data
                    del cur_rhog
                del sub_gc, rhog_ref

                if gc_user is not None:
                    del T_USER
                else:
                    del T_SWIR1, T_SWIR2, use_swir1
                del Rf_sen, omega, muv, mus

    ## add rhog to auto_grouping and update gatts
    if setu['glint_write_rhog_all'] & ('rhog' not in gem.gatts['auto_grouping']):
        gem.gatts['auto_grouping'] += ':rhog:rhogu'
    gem.gatts['glint_correction'] = 'default'
    if write: gem.gatts_update()

    ## close and return
    if (opened) & (write):
        gem.close()
        del gem
        return(gemf)
    else:
        return
