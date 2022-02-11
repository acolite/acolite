## run ACSTAR3 DSF atmospheric correction with adjacency correction
## version with minimisation of RMSD between corrected rdark and romix
## QV 2021-07-10
##                2021-07-15 (QV) added sensor specific raster psf, added force_model and force_aot
##                2021-10-06 (QV) new version with fitting and iterative methods, new psf raster/fit loading
##                2021-10-11 (QV) remove L8 PAN band from fit consideration
##                2021-10-12 (QV) added "target" fit option, added interface reflectance option
##                2021-10-14 (QV) new version of ACSTAR3
##                2021-10-21 (QV) added PSF extent check and masking of scene edges where PSF coverage is not complete

def acstar3(ncf, output=None, settings=None,
            ex = 3, ## radius of psf
            psf_raster = True,
            tolerance = 10e-4, ## tolerance of rmsd fit
            select_model = False,
            force_model = None,
            force_aot = None,
            max_wavelength = 720,
            max_iter = 20,
            luts = None,
            setu = None,
            write_rhosu = True, write_rhoa = True, write_rhoe = True,
            fit_all_bands = True, mask_edges = True,
            method = 'min_rmsd',
            target_rhos = None, target_lat = None, target_lon = None,
            include_interface = False, interface_wind = 20, interface_par = 'rsky_s',
            verbosity=0):

    import acolite as ac
    import numpy as np
    from scipy import optimize
    import sys, os, time

    import matplotlib.pyplot as plt ## for diagnostic plots alone

    ## function to retrieve adjacency corrected results
    ## for optimizing aot based on corrected rho_dark:rho_path correspondence
    ## QV 2021-07-08
    ## edited 2021-10-06 QV removed raster loading from this function
    def fit_aot_adj(aot, lut, return_results = False, return_negatives = False, plot_results = False):
        if include_interface:
            rgi_pos = [pressure, 0, raa, vza, sza, interface_wind, aot]
            rgi_pos[1] = lutdw[lut]['ipd'][interface_par]
            int_cur = {b: float(lutdw[lut]['rgi'][b](rgi_pos)) for b in bands}
        else:
            rgi_pos = [pressure, 0, raa, vza, sza, aot]


        ## compute atmosphere parameters
        atm = {'tt_gas': tg_dict['tt_gas']}
        for par in lut_par:
            rgi_pos[1] = lutdw[lut]['ipd'][par]
            atm[par] = {b: float(lutdw[lut]['rgi'][b](rgi_pos)) for b in bands}

        ## get gas optical depth from transmittance
        atm['tgas'] = {b: -np.log(atm['tt_gas'][b]) for b in bands}

        ## compute upward direct and diffuse total transmittances
        atm['udirt'] = {b: np.exp(-(atm['ttot'][b]+atm['tgas'][b])/cos_vza) for b in bands}
        atm['udift'] = {b: atm['utott'][b] - atm['udirt'][b] for b in bands}

        ## compute aerosol and rayleigh upward direct and diffuse transmittances
        atm['udifa'] = {b: atm['utott'][b] * np.exp(-0.5 * atm['tray'][b] / cos_vza) - atm['udirt'][b] for b in bands}
        atm['udifr'] = {b: atm['utott'][b] - atm['udirt'][b] - atm['udifa'][b] for b in bands}

        ## compute estimate of rhos
        data, data_mean = {}, {}
        for b in datat:
            data[b] = datat[b]/atm['tt_gas'][b] - atm['romix'][b]
            data[b] = (data[b]) / (atm['dtott'][b]*atm['utott'][b] + atm['astot'][b]*data[b])
            data_mean[b] = np.nanmean(data[b])

        ## compute psf
        psf_a, psf_cv = {}, {}
        for b in datat:
            wave = bands[b]['wave_name']
            ## raster psf
            if psf_raster_sub:
                ## load aerosol and Rayleigh nPSF from preloaded data
                ai, aw = ac.shared.closest_idx(apsfs_rasters[lut]['awv'], wave)
                pa, aa = apsfs_rasters[lut]['data']['ndfpsf_{}'.format(aw)]
                pr, ra = apsfs_rasters['Rayleigh']['data']['ndfpsf_1013'] ## use normal pressure for now

                ## crop pa and pr
                cdim = int((ex/(resolution/2000)-1)/2)
                mid = int(pa.shape[0]/2)
                pa_crop = pa[mid-cdim:mid+cdim+1, mid-cdim:mid+cdim+1]
                pr_crop = pr[mid-cdim:mid+cdim+1, mid-cdim:mid+cdim+1]

                ## new 2021-10-14
                wm = (atm['udifa'][b] * pa_crop + atm['udifr'][b] * pr_crop) / (atm['udifa'][b]+atm['udifr'][b])

                ## commented 2021-10-14
                ##
                #nPSF_diff = (atm['udifa'][b] * pa_crop + atm['udifr'][b] * pr_crop) / (atm['udifa'][b]+atm['udifr'][b])
                #nPSF_diff_cv = nPSF_diff.sum()

                #nPSF_diff *= (atm['udift'][b] * atm['dtott'][b] + data_mean[b] * atm['astot'][b]) + \
                #            ((1 - nPSF_diff_cv) * atm['udift'][b] * atm['dtott'][b] / nPSF_diff.shape[0]**2)
                ##idp = int(((nPSF_diff.shape[0]) - 1) / 2)
                #nPSF_diff[idp, idp] += atm['udirt'][b] * atm['dtott'][b]
                ##psf_a[b] = nPSF_diff

                #psf_cv[b] = wm.sum()
                #psf_a[b] = wm
                #psf_a[b] *= atm['udift'][b] * atm['dtott'][b]
                #psf_a[b][idp, idp] += atm['udirt'][b] * atm['dtott'][b]
                #psf_a[b] /= (1 - data_mean[b] * atm['astot'][b])
            ## fit psf
            else:
                coef_ray = np.array([apsfs_fits['Rayleigh'][c] for c in ['c1', 'c2', 'c3', 'c4', 'c5', 'c6']])
                coef_aer = np.array([apsfs_fits[lut][b][c] for c in ['c1', 'c2', 'c3', 'c4', 'c5', 'c6']])
                wm, exm = ac.adjacency.acstar3.w_kernel(coef_aer, coef_ray, atm['udifr'][b], atm['udifa'][b],
                                                        ex = ex, res = resolution/1000, pressure=pressure)
                ## commented 2021-10-14
                ##idp = int(((wm.shape[0]) - 1) / 2)
                #psf_cv[b] = wm.sum()
                #psf_a[b] = wm

                ## from AC v0.9 code
                # psf_a[b] = wm
                #psf_a[b] *= (atm['udift'][b] * atm['dtott'][b] + data_mean[b] * atm['astot'][b]) + \
                #            ((1 - psf_cv[b]) * atm['udift'][b] * atm['dtott'][b] / wm.shape[0]**2)
                #idp = int(((wm.shape[0]) - 1) / 2)
                #psf_a[b][idp, idp] += atm['udirt'][b] * atm['dtott'][b]

                ## from AC v1.0 code 2021-10-14
                #psf_a[b] *= (atm['udift'][b] * atm['dtott'][b] + data_mean[b] * atm['astot'][b]) + \
                #            ((1 - psf_cv[b]) * atm['udift'][b] * atm['dtott'][b] / wm.shape[0]**2)
                ###idp = int(((wm.shape[0]) - 1) / 2)
                #psf_a[b][idp, idp] += atm['udirt'][b] * atm['dtott'][b]

                #psf_a[b] *= atm['udift'][b] * atm['dtott'][b]
                #psf_a[b][idp, idp] += atm['udirt'][b] * atm['dtott'][b]
                #psf_a[b] /= (1 - data_mean[b] * atm['astot'][b])

            ## check image extent
            if (wm.shape[0] > data[b].shape[0]) | (wm.shape[1] > data[b].shape[1]):
                print('Image size {}x{} too small for PSF extent {}x{}'.format(data[b].shape[0],data[b].shape[1], wm.shape[0], wm.shape[1]))
                print('Exiting adjacency correction.')
                return

            ## new 2021-10-14
            idp = int(((wm.shape[0]) - 1) / 2)
            psf_cv[b] = wm.sum()
            psf_a[b] = wm

            psf_a[b] *= atm['udift'][b] * atm['dtott'][b]
            psf_a[b][idp, idp] += atm['udirt'][b] * atm['dtott'][b]
            psf_a[b] /= (1 - data_mean[b] * atm['astot'][b])

        ## make psf same size
        if False:
            ## get size of largest psf
            pdim = 0,0
            for b in data:
                if(psf_a[b].shape) > pdim: pdim = psf_a[b].shape
            idp_ = int(((pdim[0]) - 1) / 2)
            ## make same sizes for all psf
            for b in datat:
                cur = np.zeros(pdim)
                os = psf_a[b].shape
                hdiff = int((pdim[0]-os[0])/2), int((pdim[1]-os[1])/2)
                cur[hdiff[0]:hdiff[0]+os[0], hdiff[1]:hdiff[1]+os[1]] = psf_a[b]
                psf_a[b] = cur

        ## compute otf
        otf_a = {}
        ## make otf
        for b in datat:
            otf_a[b] = ac.adjacency.acstar3.psf_otf(data[b].shape, psf_a[b])

        ## checked with R code
        datac, datac_mean = {}, {}
        datati, datati_mean = {}, {}
        max_neg = 0

        for b in data:
            ## ipd_ varied per band
            idp_ = int(((psf_a[b].shape[0]) - 1) / 2)
            #tmp = datat[b]/atm['tt_gas'][b]
            #tmp[np.isnan(tmp)] = 0
            #ext = ac.adjacency.acstar3.extend(tmp, idp_)
            mask = np.isnan(datat[b])

            if include_interface & (interface_par == 'rsky_t'):
                ext = ac.adjacency.acstar3.extend((datat[b]/atm['tt_gas'][b])-int_cur[b], idp_, fill_nan=True)
            else:
                ext = ac.adjacency.acstar3.extend(datat[b]/atm['tt_gas'][b], idp_, fill_nan=True)

            dim_ext = ext.shape
            dsize = dim_ext[0] * dim_ext[1]

            ## fft the image data
            x_f = np.fft.fftn(ext)

            ## subtract path reflectance from first element of fft
            x_f[0,0] = x_f[0,0] - atm['romix'][b] * dsize

            #x_f[1, 1, band] <- x_f[1, 1, band] - (x_eff[band] * (1 - psf_cv[band]) *
            #   udift[band] * dtott[band] / (1 - x_eff[band] * astot[band])) * dsize
            ## new 2021-10-14
            x_f[1,1] -= (data_mean[b] * (1 - psf_cv[b]) * \
                         atm['udift'][b] * atm['dtott'][b] / (1 - data_mean[b] * atm['astot'][b])) * dsize

            ## if otf not yet transformed
            #tmp = x_f / np.fft.fftn(otf_a[b])

            ## otf is transformed in acstar3.psf_otf
            tmp = x_f / otf_a[b]
            tmp2 = np.fft.ifftn(tmp)
            tmp3 = (np.sign(tmp2) * np.abs(tmp2)).real

            ## subset extended data back to image extent
            x_s = tmp3[idp_:dim_ext[0]-idp_, idp_:dim_ext[1]-idp_]
            datac[b] = x_s
            if include_interface & (interface_par == 'rsky_s'):
                datac[b] -= int_cur[b]

            ## mask and compute mean
            datac[b][mask] = np.nan
            datac_mean[b] = np.nanmean(datac[b])

            ## mask the edges of the scene not fully covered by the PSF
            if mask_edges:
                tmp = datac[b] * np.nan
                idp_h = int(idp_/2)
                tmp[idp_h:-1*idp_h, idp_h:-1*idp_h] = datac[b][idp_h:-1*idp_h, idp_h:-1*idp_h]
                datac[b] = tmp
                del tmp

            neg = len(np.where(datac[b]<0)[0])/(datac[b].shape[0]*datac[b].shape[1])
            max_neg = max(max_neg, neg)
            ## convert back to toa w/o adjacency
            datati[b] = atm['romix'][b] + ((datac[b] * atm['udirt'][b] * atm['dtott'][b]) /
                                           (1 - datac_mean[b] * atm['astot'][b]))

        ## get rhod per band
        rhodi = {}
        for b in datat:
            if setu['dsf_spectrum_option'] == 'darkest':
                rd = np.array((np.nanpercentile(datati[b], 0)))
            if setu['dsf_spectrum_option'] == 'percentile':
                rd = np.array((np.nanpercentile(datati[b], setu['dsf_percentile'])))
            if setu['dsf_spectrum_option'] == 'intercept':
                rd = ac.shared.intercept(datati[b], setu['dsf_intercept_pixels'])
            rhodi[b] = rd

        ## compute rmsd
        rhodi_arr = np.asarray([rhodi[b] for b in rhodi])
        romix_arr = np.asarray([atm['romix'][b] for b in rhodi])
        bands_wave = [bands[b]['wave_nm'] for b in rhodi]


        if fit_all_bands:
            rmsd_all = np.sqrt(np.nanmean(np.square((rhodi_arr-romix_arr))))
        else:
            diff_sorted = np.argsort(rhodi_arr-romix_arr)[0:3]
            rmsd_all = np.sqrt(np.nanmean(np.square((rhodi_arr[diff_sorted]-romix_arr[diff_sorted]))))
            #print(diff_sorted, np.asarray([b for b in rhodi])[diff_sorted])

        ## compute rmsd with target pixel
        if method == 'target':
            target_bands = [b for b in target_rhos if b in datac]
            tar = np.asarray([target_rhos[b] for b in target_bands])
            ## was used for rsky_s
            #if add_interface:
            #    res = np.asarray([datac[b][pixi, pixj]-int_cur[b] for b in target_bands])
            #else:
            #    res = np.asarray([datac[b][pixi, pixj] for b in target_bands])
            res = np.asarray([datac[b][pixi, pixj] for b in target_bands])
            rmsd_all = np.sqrt(np.nanmean(np.square((tar-res))))

        if plot_results:
            rhod_arr = np.asarray([rhod[b] for b in rhodi])
            idx = np.argsort(bands_wave)
            oplot = ofile.replace('.nc', '_fit_results.png')
            fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
            #plt.sca(ax[0])
            plt.plot([bands_wave[i] for i in idx], [rhod_arr[i] for i in idx], '.--', color='Grey', label=r'$\rho_{dark}$ original')
            plt.plot([bands_wave[i] for i in idx], [rhodi_arr[i] for i in idx], '.--', color='Black', label=r'$\rho_{dark}$ corrected')
            plt.plot([bands_wave[i] for i in idx], [romix_arr[i] for i in idx], '.-', color='Blue', label=r'$\rho_{path}$')
            plt.xlabel('Wavelength (nm)')
            plt.ylabel(r'$\rho$ (-)')
            subt = '\n Method "{}" PSF-{}'.format(method, 'Raster' if psf_raster_sub else 'Fit')
            plt.title(r'{} $\tau_a$={:.4f} RMSD={:.2e} <0={:.1f}% {}'.format(lut[-4:], float(aot), rmsd_all, max_neg*100, subt))
            plt.legend()
            plt.savefig(oplot, dpi=300, bbox_inches='tight')
            plt.close()

            oplot = ofile.replace('.nc', '_fit_results_diff.png')
            fig, ax = plt.subplots(1,1, sharex=True, sharey=True)
            #plt.sca(ax[0])
            plt.plot([bands_wave[i] for i in idx], [rhod_arr[i]-romix_arr[i] for i in idx], '.--', color='Grey', label=r'$\rho_{dark}$ original - $\rho_{path}$')
            plt.plot([bands_wave[i] for i in idx], [rhodi_arr[i]-romix_arr[i] for i in idx], '.--', color='Black', label=r'$\rho_{dark}$ corrected - $\rho_{path}$')
            xlim = plt.xlim()
            plt.plot(xlim, (0,0), ':', color='Black')
            plt.xlim(xlim)
            plt.xlabel('Wavelength (nm)')
            plt.ylabel(r'$\rho$ difference (-)')
            subt = '\n Method "{}" PSF-{}'.format(method, 'Raster' if psf_raster_sub else 'Fit')
            plt.title(r'{} $\tau_a$={:.4f} RMSD={:.2e} <0={:.1f}% {}'.format(lut[-4:], float(aot), rmsd_all, max_neg*100, subt))
            plt.legend()
            plt.savefig(oplot, dpi=300, bbox_inches='tight')
            plt.close()


        print(lut, 'tau_a {:.5f}'.format(float(aot)), 'rmsd {:.5f}'.format(rmsd_all), 'neg {:.1f}%'.format(max_neg*100))
        if not return_results:
            if return_negatives:
                return(rmsd_all, max_neg)
            else:
                return(rmsd_all)
        else:
            return(atm, data, datac, datati, idp_)

    ## read global attributes
    gatts = ac.shared.nc_gatts(ncf)
    nc_projection = ac.shared.nc_read_projection(ncf)
    sensor = gatts['sensor']
    if (sensor not in ['L8_OLI', 'S2A_MSI', 'S2B_MSI']) and ('PlanetScope' not in sensor):
        print('ACSTAR3 not implemented for {}'.format(sensor))
        return([ncf])

    ## placeholder settings
    ## if called directly without settings
    if setu is None:
        setu = {}
        setu['dsf_spectrum_option'] = 'intercept'
        setu['dsf_spectrum_option'] = 'percentile'
        setu['dsf_percentile'] = 5
        setu['dsf_intercept_pixels'] = 200
        setu['min_tgas_aot'] = 0.85

        setu['ancillary_data'] = True
        setu['uoz_default'] = 0.3
        setu['uwv_default'] = 1.5
        setu['pressure'] = 1013.25

    ## get output from settings
    if 'output' in setu: output = setu['output']

    ## set user settings
    if 'acstar3_method' in setu: method = setu['acstar3_method']
    if 'acstar3_psf_raster' in setu: psf_raster = setu['acstar3_psf_raster']
    if 'acstar3_max_wavelength' in setu: max_wavelength = setu['acstar3_max_wavelength']
    if 'acstar3_fit_all_bands' in setu: fit_all_bands = setu['acstar3_fit_all_bands']
    if 'acstar3_write_rhosu' in setu: write_rhosu = setu['acstar3_write_rhosu']
    if 'acstar3_write_rhoa' in setu: write_rhoa = setu['acstar3_write_rhoa']
    if 'acstar3_write_rhoe' in setu: write_rhoe = setu['acstar3_write_rhoe']
    if 'acstar3_ex' in setu:
        if setu['acstar3_ex'] is not None:
            ex = float(setu['acstar3_ex'])
    if 'acstar3_mask_edges' in setu: mask_edges = setu['acstar3_mask_edges']

    ## check method
    if (method == 'target') & ((target_lat is None) | (target_lon is None) | (target_rhos is None)):
        print('Fitting method "target" requires target_lat (float), target_lon (float) and target_rhos (dict with rhos per band) to be set')
        return

    if method not in ['min_rmsd', 'iter', 'target']:
        print('Fitting method "{}" not configured, using default "iter"'.format(method))
        method = 'iter'

    ## path reflectance parameter
    rpar = 'romix'

    ## store original results
    ac_aot_550_dsf = gatts['ac_aot_550'] if 'ac_aot_550' in gatts else np.nan
    ac_model_dsf = gatts['ac_model'] if 'ac_model' in gatts else 'None'
    ac_rmsd_dsf = gatts['ac_rmsd'] if 'ac_rmsd' in gatts else np.nan

    ## get datasets and sensor rsr
    datasets = ac.shared.nc_datasets(ncf)
    ss = sensor.lower().split('_')

    resolution = gatts['scene_pixel_size'][0]
    rsrd = ac.shared.rsr_dict(sensor=sensor)

    ## get geometry
    sza, vza = gatts['sza'], gatts['vza']
    raa = gatts['raa']
    cos_vza = np.cos(vza*np.pi/180)

    ## anc data?
    if setu['ancillary_data']:
        clon = np.nanmedian(ac.shared.nc_data(ncf, 'lon'))
        clat = np.nanmedian(ac.shared.nc_data(ncf, 'lat'))
        anc = ac.ac.ancillary.get(gatts['isodate'], clon, clat)

    if 'uoz' not in gatts:
        uoz = setu['uoz_default'] * 1.0
        if setu['ancillary_data']: uoz = anc['ozone']['interp']/1000. ## convert from MET data
    else: uoz = gatts['uoz']
    if 'uwv' not in gatts:
        uwv = setu['uwv_default'] * 1.0
        if setu['ancillary_data']: uwv = anc['p_water']['interp']/10. ## convert from MET data
    else: uwv = gatts['uwv']
    if 'pressure' not in gatts:
        pressure = setu['pressure'] * 1.0
        if setu['ancillary_data']: pressure = anc['press']['interp']
    else: pressure = gatts['pressure']

    ## find target pixel
    if method == 'target':
        lat = ac.shared.nc_data(ncf, 'lat')
        lon = ac.shared.nc_data(ncf, 'lon')
        tmp = ((lon-target_lon)**2 + (lat-target_lat)**2)**0.5
        pixi, pixj = np.where(tmp == np.nanmin(tmp))
        pixi = pixi[0]
        pixj = pixj[0]
        #print(lat[pixi, pixj], lon[pixi, pixj])
        lat = None
        lon = None

    ## read LUT with required parameters
    #add_interface = include_interface & (method == 'target')
    lut_par = ['romix', 'taer', 'tray', 'utott', 'dtott', 'astot', 'ttot']
    lutdw = ac.aerlut.import_luts(add_rsky=include_interface, add_dutott=False,
                                  sensor=sensor, lut_par = lut_par) # base_luts = [lut],
    luts = list(lutdw.keys())

    ## get gas transmittance
    tg_dict = ac.ac.gas_transmittance(sza, vza, uoz=uoz, uwv=uwv, rsr=rsrd[sensor]['rsr'])

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

    ## read all apsfs fits
    apsfs_fits = {}
    for lut in luts:
        if lut[-1] == '1': model_name = 'continental'
        elif lut[-1] == '2': model_name = 'maritime'
        elif lut[-1] == '3': model_name = 'urban'
        else: continue

        ## read fit data
        rdata, adata = ac.adjacency.acstar3.read_apsfs_fits(model_name, sensor)
        if 'Rayleigh' not in apsfs_fits: apsfs_fits['Rayleigh'] = rdata
        apsfs_fits[lut] = adata

        ## temporary fix
        if len(adata) == 0:
            ref_sen = 'S2A_MSI'
            print('Warning: {} substituting closest bands from {} for {}'.format(lut, ref_sen, sensor))
            rdata, adata = ac.adjacency.acstar3.read_apsfs_fits(model_name, ref_sen)
            ##
            rsrd_ref  = ac.shared.rsr_dict(sensor=ref_sen)
            adata_e = {}
            for b in rsrd[sensor]['rsr_bands']:
                bc, bd = None, 1.
                for br in rsrd_ref[ref_sen]['rsr_bands']:
                    cbd =  rsrd_ref[ref_sen]['wave_mu'][br] - rsrd[sensor]['wave_mu'][b]
                    if np.abs(cbd) < np.abs(bd):
                        bd = cbd * 1.0
                        bc = br
                print(sensor, b, rsrd[sensor]['wave_mu'][b])
                print(ref_sen, bc, rsrd_ref[ref_sen]['wave_mu'][bc])
                adata_e[b] = adata[br]
            apsfs_fits[lut] = adata_e

    ## load psf raster files
    if psf_raster:
        apsfs_rasters = {}
        for lut in luts:
            if lut[-1] == '1':
                m = 'C'
            if lut[-1] == '2':
                m = 'M'
            if lut[-1] == '3':
                m = 'U'

            if sensor in ['L8_OLI', 'S2A_MSI', 'S2B_MSI']:
                if sensor in ['L8_OLI']:
                    psfr_a = '{}/Shared/psf_raster/{}_{}_res030_v00_p1013_grid_{}.nc'.format(ac.config['data_dir'], ss[1], ss[0], m)
                    psfr_r = '{}/Shared/psf_raster/{}_{}_res030_v00_p1100-500_grid_R.nc'.format(ac.config['data_dir'], ss[1], ss[0], m)
                if sensor in ['S2A_MSI', 'S2B_MSI']:
                    psfr_a = '{}/Shared/psf_raster/{}_{}_res010_v00_p1013_grid_{}.nc'.format(ac.config['data_dir'], ss[1], ss[0], m)
                    psfr_r = '{}/Shared/psf_raster/{}_{}_res010_v00_p1100-500_grid_R.nc'.format(ac.config['data_dir'], ss[1], ss[0], m)
                ads = ac.shared.nc_datasets(psfr_a)
                awv = [ds.split('_')[1] for ds in ads if 'ndfpsf_' in ds]
                print('Loaded raster psf for {} {}'.format(sensor, lut))
                print('Aerosol: {}'.format(psfr_a))
                apsfs_rasters[lut] = {'ads': ads, 'awv': awv}
                apsfs_rasters[lut]['data'] = {ds: ac.shared.nc_data(psfr_a, ds, attributes=True) for ds in ads if 'ndfpsf_' in ds}
                if 'Rayleigh' not in apsfs_rasters:
                    print('Rayleigh: {}'.format(psfr_r))
                    rds = ac.shared.nc_datasets(psfr_r)
                    apsfs_rasters['Rayleigh'] = {'data': {ds: ac.shared.nc_data(psfr_r, ds, attributes=True) for ds in rds if 'ndfpsf_' in ds}}
            else:
                print('Sensor raster file not available {}'.format(sensor))
        if len(apsfs_rasters) == 0:
            psf_raster = False

    psf_raster_sub = psf_raster

    ## read toa data - for bands used in the aot optimisation
    datat = {}
    for b in bands:
        if tg_dict['tt_gas'][b] < setu['min_tgas_aot']: continue
        if bands[b]['rhot_ds'] not in datasets: continue
        if bands[b]['wave_nm'] > max_wavelength: continue
        if (sensor == 'L8_OLI') & (b == '8'): continue ## skip PAN band for L8, doesn't work nicely with the iterative method
        datat[b] = ac.shared.nc_data(ncf, bands[b]['rhot_ds'])
        datat[b][datat[b].mask] = np.nan
        datat[b] = datat[b].data

    ## get rhod per band - not used, could be initial estimate of aot
    if False:
        rhod = {}
        for b in datat:
            if setu['dsf_spectrum_option'] == 'darkest':
                rd = np.array((np.nanpercentile(datat[b], 0)))
            if setu['dsf_spectrum_option'] == 'percentile':
                rd = np.array((np.nanpercentile(datat[b], setu['dsf_percentile'])))
            if setu['dsf_spectrum_option'] == 'intercept':
                rd = ac.shared.intercept(datat[b], setu['dsf_intercept_pixels'])
            rhod[b] = rd
        print(rhod)

        ## fit aod per band and model
        rhod_a = {}
        rhod_arr = np.array([rhod[b] for b in rhod])
        for li, lutn in enumerate(luts):
            rhod_a[lutn] = {'fit':{}}
            ## fit aot to each band
            for b in rhod:
                tmp = lutdw[lutn]['rgi'][b]((pressure, lutdw[lutn]['ipd'][rpar],raa,vza,sza,lutdw[lutn]['meta']['tau']))
                rhod_a[lutn]['fit'][b] = np.interp(rhod[b], tmp,lutdw[lutn]['meta']['tau'], left=np.nan, right=np.nan)

            ## select min aot and compute romix
            rhod_a[lutn]['sorted'] = sorted(rhod_a[lutn]['fit'].items(), key =lambda kv:(kv[1], kv[0]))
            rhod_a[lutn]['min'], rhod_a[lutn]['aot'] = rhod_a[lutn]['sorted'][0]

            rhod_a[lutn]['romix'] = {b: float(lutdw[lutn]['rgi'][b]((pressure, lutdw[lutn]['ipd'][rpar],raa,vza,sza,rhod_a[lutn]['aot']))) for b in rhod }
            rhod_a[lutn]['romix_arr'] = [rhod_a[lutn]['romix'][b] for b in rhod]

            ## compute rmsd
            rhod_a[lutn]['rmsd_all'] = np.sqrt(np.nanmean(np.square((rhod_arr-rhod_a[lutn]['romix_arr']))))
            rhod_a[lutn]['rmsd'] = np.sqrt(np.nanmean(np.square((np.asarray([rhod_a[lutn]['romix'][s[0]] for s in rhod_a[lutn]['sorted'][0:2]])-
                                                                np.asarray([rhod[s[0]] for s in rhod_a[lutn]['sorted'][0:2]])))))

            ## print
            print(lutn, rhod_a[lutn]['sorted'])

    ## make output directory here so we can safely output plots
    if not os.path.exists(output): os.makedirs(output)

    ## get results for both luts
    res = {}
    for lut in luts:
        if force_model is not None:
            if lut != force_model: continue
            if force_aot is not None:
                res[lut] = {'rmsd': -1, 'aot': float(force_aot)}
                continue

        ## method "fit" minimize rmsd
        if method == 'min_rmsd':
            print('Fitting {}'.format(lut))
            t0 = time.time()
            try:
                xro = optimize.minimize(fit_aot_adj, [0.001], args=(lut), bounds = [(0.001,5.0)], tol=tolerance, options={'maxiter': max_iter})
                res[lut] = {'rmsd': xro.fun, 'aot': xro.x[0], 'fit': xro}
                print('Fitting {} took {:.1f} seconds, aot550 = {:.4f}'.format(lut, time.time()-t0, xro.x[0]))
            except:
                print('Fitting failed.')
                return

        ## method "iter" iterate until conditions are met
        elif method == 'iter': ## iterative option
            print('Fitting {}'.format(lut))
            t0 = time.time()

            aots = {}

            ## initial range
            nsamples = 3
            arange = 0, lutdw[lut]['meta']['tau'][-2]
            step = (arange[1]-arange[0])/(nsamples+1)
            cb = np.linspace(arange[0], arange[1], num=nsamples+2)
            cb[0] = 0.001 ## set to 0.001

            iterate = True
            iteration = 0
            #max_iterations = 50

            ## tuneable parameters
            min_neg = 0.01 ## 1% negatives in VNIR
            min_neg = 0.005 ## 0.5% negatives in VNIR
            min_neg = 0.001 ## 0.1% negatives in VNIR
            min_step = 0.001

            while iterate:
                iteration+=1
                print(lut, iteration, arange, step)

                ## get results for current aot steps
                for aot in cb:
                    if aot in aots: continue
                    #rmsd, neg = fit_aot_adj(aot, lut, return_results = False, return_negatives = True)
                    ret = fit_aot_adj(aot, lut, return_results = False, return_negatives = True)
                    if ret is None:
                        return
                    else:
                        rmsd, neg = ret
                    aots[aot] = {'neg':neg, 'rmsd':rmsd}

                ## find aot position with minimal negatives (but we want a small amount of negatives)
                aot_list = sorted(list(aots.keys()))
                neg_list = [aots[aot]['neg'] for aot in aot_list]
                i1 = np.where(np.asarray(neg_list) >= min_neg)[0]
                if len(i1) == 0:
                    print('Found no negatives in aot range {}-{}, image is likely fully cloudy.'.format(arange[0], arange[1]))
                    iterate = False
                else:
                    i1 = i1[0]
                    i0 = np.max((0, i1-1))

                ## stop iterating
                if np.abs(step) <= min_step:
                    print('Reached min step {} <= {}'.format(step, min_step))
                    iterate = False
                if iteration>= max_iter: iterate = False

                if iterate:
                    arange = aot_list[i0], aot_list[i1]
                    step = (arange[1]-arange[0])/(nsamples+1)
                    cb = np.linspace(arange[0], arange[1], num=nsamples+2)

            aot_list = np.asarray(sorted(list(aots.keys())))
            neg_list = np.asarray([aots[aot]['neg'] for aot in aot_list])
            rmsd_list = np.asarray([aots[aot]['rmsd'] for aot in aot_list])
            i = np.argsort(abs(aot_list-cb[2]))[0]
            res[lut] = {'rmsd': rmsd_list[i], 'aot': aot_list[i], 'fit': {'fun': rmsd_list[i]}}
            print('Fitting {} took {:.1f} seconds, aot550 = {:.4f}'.format(lut, time.time()-t0, res[lut]['aot']))

        ## method "target" minimize rmsd with target pixel
        elif method == 'target':
            print('Fitting {}'.format(lut))
            t0 = time.time()
            #if add_interface:
            #    xro = optimize.minimize(fit_aot_adj, [0.001, 0], args=(lut), bounds = [(0.001,5.0), (1,20.0)], tol=tolerance)
            #    res[lut] = {'rmsd': xro.fun, 'aot': xro.x[0], 'fit': xro, 'intf': xro.x[1]}
            #    print('Fitting {} took {:.1f} seconds, aot550 = {:.4f}, intf = {:.4f}'.format(lut, time.time()-t0, xro.x[0], xro.x[1]))
            #else:
            #    xro = optimize.minimize(fit_aot_adj, [0.001], args=(lut), bounds = [(0.001,5.0)], tol=tolerance)
            #    res[lut] = {'rmsd': xro.fun, 'aot': xro.x[0], 'fit': xro}
            #    print('Fitting {} took {:.1f} seconds, aot550 = {:.4f}'.format(lut, time.time()-t0, xro.x[0]))
            try:
                xro = optimize.minimize(fit_aot_adj, [0.001], args=(lut), bounds = [(0.001,5.0)], tol=tolerance, options={'maxiter': max_iter})
                res[lut] = {'rmsd': xro.fun, 'aot': xro.x[0], 'fit': xro}
                print('Fitting {} took {:.1f} seconds, aot550 = {:.4f}'.format(lut, time.time()-t0, xro.x[0]))
            except:
                print('Fitting failed.')
                return
    #return(bands, res)

    ## read additional toa data
    for b in bands:
        if b in datat: continue
        if bands[b]['rhot_ds'] not in datasets: continue
        if tg_dict['tt_gas'][b] < setu['min_tgas_aot']: continue
        datat[b] = ac.shared.nc_data(ncf, bands[b]['rhot_ds'])
        datat[b][datat[b].mask] = np.nan
        datat[b] = datat[b].data

    ## get rhod per band - to plot original rhod
    if True:
        rhod = {}
        for b in datat:
            if setu['dsf_spectrum_option'] == 'darkest':
                rd = np.array((np.nanpercentile(datat[b], 0)))
            if setu['dsf_spectrum_option'] == 'percentile':
                rd = np.array((np.nanpercentile(datat[b], setu['dsf_percentile'])))
            if setu['dsf_spectrum_option'] == 'intercept':
                rd = ac.shared.intercept(datat[b], setu['dsf_intercept_pixels'])
            rhod[b] = rd

    ## make outputfile
    if 'rhosu' not in gatts['auto_grouping']:
        gatts['auto_grouping']+=':rhosu'
    if 'rhoe' not in gatts['auto_grouping']:
        gatts['auto_grouping']+=':rhoe'
    if 'rhoa' not in gatts['auto_grouping']:
        gatts['auto_grouping']+=':rhoa'

    ## which luts to output
    luts_out = [luts[np.argsort([res[lut]['rmsd'] for lut in res])[0]]] if select_model else list(res.keys())
    print('Outputting results for {}'.format(','.join(luts_out)))

    ## track generated files
    ofiles = []

    ## per lut
    for lut in luts_out:
        print('Writing results for {}'.format(lut))
        ## make outputfile
        if output is None:
            odir = os.path.dirname(gatts['ofile'])
        else:
            odir = '{}'.format(output)
        ofile = '{}/{}_L2R_ACSTAR3{}.nc'.format(odir, gatts['oname'], '' if len(luts_out) == 1 else '_{}'.format(lut[-4:]))

        ## retrieve atm and results
        #if 'intf' in res[lut]:
        #    atm, data, datac, datati, idp = fit_aot_adj((res[lut]['aot'],res[lut]['intf']), lut, return_results = True, plot_results = True)
        #else:
        #    atm, data, datac, datati, idp = fit_aot_adj(res[lut]['aot'], lut, return_results = True, plot_results = True)

        atm, data, datac, datati, idp = fit_aot_adj(res[lut]['aot'], lut, return_results = True, plot_results = True)

        ## store original results
        gatts['ac_aot_550_dsf'] = ac_aot_550_dsf
        gatts['ac_model_dsf'] = ac_model_dsf
        gatts['ac_rmsd_dsf'] = ac_rmsd_dsf

        ## set new results
        gatts['ac_aot_550'] = res[lut]['aot']
        gatts['ac_model'] = lut
        gatts['ac_rmsd'] = res[lut]['fit']['fun']

        new = True
        for ds in datasets:
            if ('rhot_' in ds) & (ds not in ['rhot_1373']) : continue
            if 'rhos_' in ds: continue
            if (nc_projection is not None) & (ds in ['x', 'y', gatts['projection_key']]): continue
            d, att = ac.shared.nc_data(ncf, ds, attributes = True)
            ac.output.nc_write(ofile, ds, d, dataset_attributes = att, attributes = gatts,
                                new = new, nc_projection=nc_projection)
            new = False

        ## compute adjacency and environment effects
        datac_mean = {}
        datata, datasa = {}, {}
        for b in bands:
            if b not in data: continue
            att = bands[b]
            for k in atm: att[k] = atm[k][b]

            ## write rhot
            ac.output.nc_write(ofile, att['rhot_ds'], datat[b], dataset_attributes = att)
            if atm['tt_gas'][b] < setu['min_tgas_aot']: continue ## dont compute surface for these

            ## get mask
            mask = np.isnan(datat[b])

            ## write rhos
            datac[b][mask] = np.nan
            ac.output.nc_write(ofile, att['rhos_ds'], datac[b], dataset_attributes = att)

            ## write rhosu
            if write_rhosu:
                data[b][mask] = np.nan
                ac.output.nc_write(ofile, att['rhos_ds'].replace('rhos_', 'rhosu_'), data[b], dataset_attributes = att)

            ## mean surface with no adjacency
            datac_mean[b] = np.nanmean(datac[b])

            ## adjacency at toa
            if write_rhoa:
                datata[b] = datat[b]/atm['tt_gas'][b] - datati[b]
                datata[b][mask] = np.nan
                ac.output.nc_write(ofile, att['rhos_ds'].replace('rhos_', 'rhoa_'), datata[b], dataset_attributes = att)

            ## adjacency at boa
            if write_rhoe:
                datasa[b] = datata[b] / (atm['udift'][b] * atm['dtott'][b] + datac_mean[b] * atm['astot'][b])
                datasa[b][mask] = np.nan
                ac.output.nc_write(ofile, att['rhos_ds'].replace('rhos_', 'rhoe_'), datasa[b], dataset_attributes = att)

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
            rsrd_o = ac.shared.rsr_dict(sensor_o)[sensor_o]
            ob = {k:rsrd_o[k]['O'] for k in ['wave_mu', 'wave_nm', 'wave_name']}
            ob['rhos_ds'] = 'rhos_{}'.format(ob['wave_name'])
            ob['wavelength'] = ob['wave_nm']
            for c in ob_cfg: ob[c] = ob_cfg[c]

            ## compute orange band (corrected)
            ob_data =  datac['8'] * float(ob_cfg['pf'])
            ob_data += datac['3'] * float(ob_cfg['gf'])
            ob_data += datac['4'] * float(ob_cfg['rf'])
            ac.output.nc_write(ofile, 'rhos_{}'.format(ob['wave_name']), ob_data, dataset_attributes=ob)

            ## compute orange band (uncorrected)
            ob_data =  data['8'] * float(ob_cfg['pf'])
            ob_data += data['3'] * float(ob_cfg['gf'])
            ob_data += data['4'] * float(ob_cfg['rf'])
            ac.output.nc_write(ofile, 'rhosu_{}'.format(ob['wave_name']), ob_data, dataset_attributes=ob)
            del ob_data
        ofiles.append(ofile)

    ## delete some stuff
    del lutdw
    del data
    del datac
    del datati
    del datat
    del datata
    del datasa

    return(ofiles)
