## def radcor
## main radcor function
##
## written by Alexandre Castagna, UGent (R version)
##            Quinten Vanhellemont, RBINS (Python version)
##
## for the RAdCor project 2023-2024
##
## modifications: 2024-03-20 (QV) converted to function
##                2024-03-25 (QV) added user defined and DSF derived aot/model (& moved tau_Ray forward)
##                                if radcor_force_model, run only using that model
##                2024-03-26 (QV) added TS-DSF/DSF centre crop using radcor_aot_estimate_centre_extent
##                2024-03-27 (QV) added mask of rhos output, added option to set TS-DSF centre pixel of PSF to 0 (tsdsf_psf_centre_zero=True)
##                2024-03-28 (QV) added option to mask to radcor_aot_estimate_centre_extent (radcor_aot_estimate_centre_mask=False)
##                                fixed using region_name in output file name
##                2024-04-02 (QV) changed handling of atmospheric correction to inner function correct_band
##                                added option to optimise aot to in situ measurement
##                2024-04-22 (QV) moved outputfile creation to after returns for incorrect settings
##                                handling of NetCDF through gem class
##                2024-04-24 (QV) added acolite file type and RAdCor version
##                                added check for forced model, split off import_fits function
##                2024-05-21 (QV) added aot optimise options RMSD/MARD
##                                added pixel mask for optimisation
##                2024-05-22 (QV) added romix+rsky_t option
##                2024-05-27 (QV) added elevation/dem pressure option
##                2024-06-05 (QV) fixed issue for ancillary_data=False
##                2024-09-24 (QV) added AC edits and renamed setting parameters
##                2024-10-14 (QV) added tsdsf_wave_range
##                2024-10-16 (QV) fixed output of rhot cirrus/wv and output of bands excluded by tsdsf_wave_range
##                2024-11-06 (QV) moved output file creation, added radcor_crop_centre
##                2024-11-20 (QV) added "adjacency corrected rhot" rhotc
##                2024-12-02 (QV) fixed crop for optimised results
##                2024-12-16 (QV) removed radcor/tsdsf_kernel_rescale and added renormalise to radcor/tsdsf_kernel_complete_method
##                2025-01-21 (QV) added radcor_write_rhotc_separate_file option
##                2025-01-24 (QV) added output_ed option
##                2025-02-03 (QV) use downward gas transmittance for output_ed
##                2025-02-04 (QV) updated settings parsing
##                2025-02-10 (QV) renamed radcor_optimise_* settings to optimise_* settings, added optimise_target_rhos_file
##                2025-03-17 (QV) use separate function for reading optim target
##                2025-03-19 (QV) added s3_product_type for MERIS/OLCI

def radcor(ncf, settings = None):
    import os, json
    import acolite as ac
    import numpy as np
    import scipy.optimize

    ### correct_band
    ## nested function to perform correction
    ## QV 2024-04-02 split off and adapted from main function
    ## if write is True then the output file is written (i.e. for final call with aot and lut determined)
    ## else the rho_s is returned
    def correct_band(b, aot, new = False, write = True, quiet = True):
        ## Load rhot data
        if not quiet: print('Loading TOA data for band {} ({})'.format(b, bands[b]['rhot_ds']))
        rho_toa, att = gem.data(bands[b]['rhot_ds'], attributes = True)
        rho_toa_mask = np.isnan(rho_toa)

        ## Gas transmittance correction
        rho_toa /= bands[b]['tt_gas']
        rho_toa_av = np.nanmean(rho_toa)

        ## Get atmospheric parameters
        if add_rsky:
            tau_ray = lutdw[luts[0]]['rgi'][b]((pressure, lutdw[luts[0]]['ipd']['tray'], raa, vza, sza, wind, 0.001))
            rho_a     = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd'][romix_par], raa, vza, sza, wind, aot))
            rho_a_sph = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['astot'], raa, vza, sza, wind, aot))
            T_d_tot   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['dtott'], raa, vza, sza, wind, aot))
            T_u_tot   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['utott'], raa, vza, sza, wind, aot))
            tau_tot   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['ttot'],  raa, vza, sza, wind, aot))
        else:
            tau_ray = lutdw[luts[0]]['rgi'][b]((pressure, lutdw[luts[0]]['ipd']['tray'], raa, vza, sza, 0.001))
            rho_a     = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['romix'], raa, vza, sza, aot))
            rho_a_sph = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['astot'], raa, vza, sza, aot))
            T_d_tot   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['dtott'], raa, vza, sza, aot))
            T_u_tot   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['utott'], raa, vza, sza, aot))
            tau_tot   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['ttot'],  raa, vza, sza, aot))
        T_u_dir   = np.exp(-tau_tot / cos_vza)
        T_u_dif   = T_u_tot - T_u_dir

        # Here we are already removing rho_a, so it should never get the previous values.
        x_a = 1.0 * rho_toa - rho_a
        ## fill nans
        x_a[np.isnan(x_a)] = rho_toa_av - rho_a

        ## Extend rhot array if requested to keep edges
        if setu['radcor_edge_extend']:
            x_a_ = ac.adjacency.radcor.extend.extend(x_a, id_psf, x_a_dim, x_f_dim, method = setu['radcor_edge_extend_method'])
            ## Do FFT of rho_toa data
            if not quiet: print('Computing fft of extended TOA array')
            x_f = np.fft.fftn(x_a_)
        else:
            ## Do FFT of rho_toa data
            if not quiet: print('Computing fft of TOA array')
            x_f = np.fft.fftn(x_a)

        ## Get PSF for selected model and current band
        aer_psf  = coefs_psf_aer[am][b]
        coef_psf_aer = [aer_psf['coefficients'][c] for c in aer_psf['coefficients']]
        coef_psf_ray = [coefs_psf_ray['fa_generic_R']['coefficients'][c] for c in coefs_psf_ray['fa_generic_R']['coefficients']]

        ## Get SAF for selected model and current band
        aer_saf  = coefs_saf_aer[am][b]
        coef_saf_aer = [aer_saf['coefficients'][c] for c in aer_saf['coefficients']]
        coef_saf_ray = [coefs_saf_ray['fa_generic_R']['coefficients'][c] for c in coefs_saf_ray['fa_generic_R']['coefficients']]

        ## Compute mixed SAF and mixed PSF
        if add_rsky:
            T_u_ray = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['utotr'], raa, vza, sza, wind, aot))
            T_u_aer = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['utota'], raa, vza, sza, wind, aot))
            tau_aer = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['taer'],  raa, vza, sza, wind, aot))
            rho_a_sph_ray = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['asray'], raa, vza, sza, wind, aot)) ## AC 2024-09-09 ## In order to properly weight the SAF
            rho_a_sph_aer = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['asaer'], raa, vza, sza, wind, aot)) ## AC 2024-09-09 ##
        else:
            T_u_ray = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['utotr'], raa, vza, sza, aot))
            T_u_aer = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['utota'], raa, vza, sza, aot))
            tau_aer = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['taer'],  raa, vza, sza, aot))
            rho_a_sph_ray = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['asray'], raa, vza, sza, aot)) ## AC 2024-09-09 ## In order to properly weight the SAF
            rho_a_sph_aer = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['asaer'], raa, vza, sza, aot)) ## AC 2024-09-09 ##

        T_u_dif_ray = T_u_ray - np.exp(-tau_ray / cos_vza)
        T_u_dif_aer = T_u_aer - np.exp(-tau_aer / cos_vza)

        psf_mix, ext = ac.adjacency.radcor.tools.w_kernel(coef_psf_aer, coef_psf_ray, T_u_dif_ray, T_u_dif_aer,
            ex = setu['radcor_kernel_radius'], res = resolution / 1000, pressure = pressure)
        #saf_mix, ext = ac.adjacency.radcor.tools.w_kernel(coef_saf_aer, coef_saf_ray, T_u_dif_ray, T_u_dif_aer,
        #    ex = setu['radcor_kernel_radius'], res = resolution / 1000, pressure = pressure) ## AC 2024-09-09 ## In order to properly weight the SAF
        saf_mix, ext = ac.adjacency.radcor.tools.w_kernel(coef_saf_aer, coef_saf_ray, rho_a_sph_ray, rho_a_sph_aer,
            ex = setu['radcor_kernel_radius'], res = resolution / 1000, pressure = pressure) ## AC 2024-09-09 ## rho_a_sph is not equal to rho_a_sph_ray + rho_a_sph_aer, but seems a better approximation than using Upward transmittance

        saf_cv_mix = np.sum(saf_mix)
        psf_cv_mix = np.sum(psf_mix)

        if (psf_cv_mix > 1):
            if not quiet: print("Integro-normalized PSF integral greater than 1.")#Reduce the PSF extent.")
        if (saf_cv_mix > 1):
            if not quiet: print("Integro-normalized SAF integral greater than 1.")#Reduce the PSF extent.")

        if not quiet: print('PSF cov: {}'.format(psf_cv_mix)) ## AC 2024-09-09 ## Moved here, to provide an idea of the actual coverage before transformations

        if (setu['radcor_kernel_complete_method'] == 'renormalise'):
            saf_mix /= saf_cv_mix
            saf_cv_mix = 1.0
            psf_mix /= psf_cv_mix
            psf_cv_mix = 1.0
        elif (setu['radcor_kernel_complete_method'] == "neighbourhood"): ## AC 2024-09-09 ## In order to avoid subtractions when using the neighbourhood method
            saf_mix += (1.0 - saf_cv_mix) * psf_uni ## AC 2024-09-09 ##
            saf_cv_mix = 1.0 ## AC 2024-09-09 ##
            psf_mix += (1.0 - psf_cv_mix) * psf_uni ## AC 2024-09-09 ##
            psf_cv_mix = 1.0 ## AC 2024-09-09 ##

        #if not quiet: print('PSF cov: {}'.format(psf_cv_mix)) ## AC 2024-09-09 ## Moved up, to provide an idea of the actual coverage before transformations

        psf_mix_sc = psf_mix * T_u_dif * T_d_tot
        psf_mix_sc[id_psf, id_psf] += (T_u_dir * T_d_tot)

        ## START DEVELOPMENT BLOCK ##
        if not quiet:
            if setu['radcor_development']:
                a_ = {'units': '1', '_FillValue': 1e+30, 'long_name': 'Discrete normalized mixed SAF',
                      'parameter': 'psf_ray_' + bands[b]['wave_name'], 'standard_name': 'normalized_mixed_saf',
                      'wavelength': bands[b]['wave_name'], 'coverage': saf_cv_mix}
                gempsf.write('saf_mix_' + bands[b]['wave_name'], saf_mix, ds_att = a_)

                saf_mix_cov[bi] = saf_cv_mix
                a_ = {'units': '1', '_FillValue': 1e+30, 'long_name': 'Discrete normalized mixed diffuse PSF',
                      'parameter': 'psf_mix_' + bands[b]['wave_name'], 'standard_name': 'normalized_mixed_diffuse_psf',
                      'wavelength': bands[b]['wave_name'], 'coverage': psf_cv_mix}
                gempsf.write('psf_mix_' + bands[b]['wave_name'], psf_mix, ds_att = a_)

                psf_mix_cov[bi] = psf_cv_mix
                a_ = {'units': '1', '_FillValue': 1e+30, 'long_name': 'Discrete scaled mixed total PSF',
                      'parameter': 'psf_scl_' + bands[b]['wave_name'], 'standard_name': 'scaled_mixed_total_psf',
                      'wavelength': bands[b]['wave_name'], 'coverage': psf_cv_mix}
                gempsf.write('psf_scl_' + bands[b]['wave_name'], psf_mix_sc, ds_att = a_)
        ## END DEVELOPMENT BLOCK ##

        ## Compute OTF
        otf_saf = np.zeros(x_f_dim)
        otf_saf[0:saf_mix.shape[0], 0:saf_mix.shape[1]] = 1.0 * saf_mix
        otf_saf = np.fft.fftn(otf_saf)

        otf_mix = np.zeros(x_f_dim)
        otf_mix[0:psf_mix.shape[0], 0:psf_mix.shape[1]] = 1.0 * psf_mix
        otf_mix = np.fft.fftn(otf_mix)

        otf_mix_sc = np.zeros(x_f_dim)
        otf_mix_sc[0:psf_mix_sc.shape[0], 0:psf_mix_sc.shape[1]] = 1.0 * psf_mix_sc
        otf_mix_sc = np.fft.fftn(otf_mix_sc)

        # All rho_env_psf calculation has been moved to the end of the band correction. ## AC 2024-09-13 ##
        #rho_env_psf = ac.adjacency.radcor.tools.conv(x_f, otf_mix, id_psf, x_f_dim, keep_edges = setu['radcor_edge_extend']) ## AC 2024-09-09 ## This is not necessary here, and it is not the correct estimate of rho_env_psf. It should be done from rho_s_est, rescaling for the spherical albedo effect after the surface estimation. However, it will be complex for me to write this section of the code and since is not relevant for the research now, I will leave it as it is.

        ## Estimate neighbourhood
        #if (not setu['radcor_kernel_rescale']) & (setu['radcor_kernel_complete_method'] == "neighbourhood"): ## AC 2024-09-09 ## No longer performed as psf_uni was already added to psf_mix
        #    rho_nbh_toa = ac.adjacency.radcor.tools.conv(x_f, otf_uni, id_psf, x_f_dim, keep_edges = setu['radcor_edge_extend']) ## AC 2024-09-09 ##

        #if (not setu['radcor_kernel_rescale']): ## AC 2024-09-09 ## Note the comment on rho_env_psf
        #    if (setu['radcor_kernel_complete_method'] == "average"):
        #        rho_env_psf += (1 - psf_cv_mix) * rho_toa_av
        #    #if (setu['radcor_kernel_complete_method'] == "neighbourhood"): ## AC 2024-09-09 ## No longer performed as psf_uni was already added to psf_mix
        #    #    rho_env_psf += (1 - psf_cv_mix) * rho_nbh_toa ## AC 2024-09-09 ##

        #rho_env_psf = rho_env_psf / (T_d_tot * T_u_tot)
        #rho_env_psf_av = np.nanmean(rho_env_psf) ## AC 2024-09-09 ## This is not necessary, as "average" complete should not use this average, but the scene average.

        ## Estimate surface reflectance with spherical albedo
        rho_env_sph_toa_est = ac.adjacency.radcor.tools.conv(x_f, otf_saf, id_psf, x_f_dim, keep_edges = setu['radcor_edge_extend'])

        if (setu['radcor_kernel_complete_method'] == "average"):
            rho_env_sph_toa_est += (1 - saf_cv_mix) * rho_toa_av
        #if (setu['radcor_kernel_complete_method'] == "neighbourhood"): ## AC 2024-09-09 ## No longer performed as psf_uni was already added to saf_mix
        #    rho_env_sph_toa_est += (1 - saf_cv_mix) * rho_nbh_toa ## AC 2024-09-09 ##

        rho_env_sph_est = rho_env_sph_toa_est / (T_d_tot * T_u_tot + rho_env_sph_toa_est * rho_a_sph)

        ## START DEVELOPMENT BLOCK ##
        if not quiet:
            if setu['radcor_development']:
                a_ = {'units': '1', '_FillValue': 1e+30, 'long_name': 'rho_env_sph estimate',
                      'parameter': 'rho_hyp_' + am + '_' + bands[b]['wave_name'], 'standard_name': 'rho_env_sph_estimate',
                      'wavelength': bands[b]['wave_name']}
                gemconv.write('rho_sph_' + am + '_' + bands[b]['wave_name'], rho_env_sph_est, ds_att = a_)

        ## END DEVELOPMENT BLOCK ##

        ## Need to recompute FFT here since we will subtract missing PSF
        if setu['radcor_kernel_complete_method'] == "average":
            #x_a -= (T_d_tot * rho_env_psf_av * (1 - psf_cv_mix) * T_u_dif) ## AC 2024-09-09 ## This should be rho_toa_av
            x_a -= (T_d_tot * rho_toa_av * (1 - psf_cv_mix) * T_u_dif) ## AC 2024-09-09 ##
        #if setu['radcor_kernel_complete_method'] == "neighbourhood": ## AC 2024-09-09 ## No longer performed as psf_uni was already added to psf_mix
            #x_a -= (rho_nbh_toa * (1 - psf_cv_mix) * T_u_dif) ## AC 2024-09-09 ##
            ## Extend rhot array if requested to keep edges
            if setu['radcor_edge_extend']: ## AC 2024-09-09 ## Moved here because it only needs to be performed for the method "average", considering that for "neighbourhood" the psf_uni was added to psf_mix and saf_mix
                x_a_ = ac.adjacency.radcor.extend.extend(x_a, id_psf, x_a_dim, x_f_dim, method = setu['radcor_edge_extend_method']) ## AC 2024-09-09 ##
                ## Do FFT of rho_toa data
                if not quiet: print('Computing fft of TOA array') ## AC 2024-09-09 ##
                x_f = np.fft.fftn(x_a_) ## AC 2024-09-09 ##
            else: ## AC 2024-09-09 ##
                ## Do FFT of rho_toa data
                if not quiet: print('Computing fft of TOA array') ## AC 2024-09-09 ##
                x_f = np.fft.fftn(x_a) ## AC 2024-09-09 ##

        ## Extend rhot array if requested to keep edges ## AC 2024-09-09 ## This block was moved inside the 'if setu['radcor_kernel_complete_method'] == "average":'
        #if setu['radcor_edge_extend']: ## AC 2024-09-09 ##
        #    x_a_ = ac.adjacency.radcor.extend.extend(x_a, id_psf, x_a_dim, x_f_dim, method = setu['radcor_edge_extend_method']) ## AC 2024-09-09 ##
        #    ## Do FFT of rho_toa data ## AC 2024-09-09 ##
        #    if not quiet: print('Computing fft of TOA array') ## AC 2024-09-09 ##
        #    x_f = np.fft.fftn(x_a_) ## AC 2024-09-09 ##
        #else: ## AC 2024-09-09 ##
        #    ## Do FFT of rho_toa data ## AC 2024-09-09 ##
        #    if not quiet: print('Computing fft of TOA array') ## AC 2024-09-09 ##
        #    x_f = np.fft.fftn(x_a) ## AC 2024-09-09 ##

        ## Deconvolution and estimate rhos
        rho_s_est = ac.adjacency.radcor.tools.deconv(x_f, otf_mix_sc, id_psf, x_f_dim, keep_edges = setu['radcor_edge_extend'])
        rho_env_psf = (rho_toa - rho_a - T_d_tot * rho_s_est * T_u_dir) / (T_d_tot * T_u_dif) ## AC 2024-09-13 ##
        rho_s_est *= (1 - rho_env_sph_est * rho_a_sph)                # AC_20240221 (it was included in the PSF scaling)

        ## add back TOA mask
        rho_s_est[rho_toa_mask] = np.nan
        rho_env_psf[rho_toa_mask] = np.nan

        if not write: ## return rhos array
            return(rho_s_est)
        else: ## write to output file
            ## mask to centre where aot was estimated  ## QV 2024-03-28
            ## indices also can correspond to full PSF coverage so do not apply if edges are extended
            if (setu['radcor_aot_estimate_centre_mask']) & (not setu['radcor_edge_extend']):
                rho_s_est[0:cen_offset_0, :] = np.nan ## top
                rho_s_est[:, 0:cen_offset_1] = np.nan ## left
                rho_s_est[x_a_dim[0] - cen_offset_0:, :] = np.nan ## bottom
                rho_s_est[:, x_a_dim[1] - cen_offset_1:] = np.nan ## right
            ## mask edges
            if not setu['radcor_edge_extend']: rho_s_est *= scene_edge_mask

            ## Compute homogeneous rhos
            if setu['radcor_write_rhosu']: # AC 20240306
                rho_s_homo = rho_toa - rho_a
                rho_s_homo /= (T_d_tot * T_u_tot + rho_s_homo * rho_a_sph)

            if not quiet: print('Writing band {}'.format(b))
            for k in bands[b]: att[k] = bands[b][k] ## add some extra

            ## write rhot
            if setu['radcor_write_rhot']:
                if setu['radcor_crop_centre']: ## crop to centre are
                    rho_toa = rho_toa[cen_offset_0:x_a_dim[0] - cen_offset_0, cen_offset_1:x_a_dim[1] - cen_offset_1]
                gemo.write(bands[b]['rhot_ds'], rho_toa, ds_att = att)

            ## write toa corrected for adjacency
            ## i.e. rhos transfered back to toa for homogeneous surface
            if setu['radcor_write_rhotc']:
                rhotc = (rho_s_est * T_d_tot * T_u_tot) / (1 - rho_a_sph * rho_s_est)
                ## add path reflectance
                rhotc +=  rho_a
                ## add back gas transmittance
                rhotc *= bands[b]['tt_gas']
                if setu['radcor_crop_centre']: ## crop to centre area
                    rhotc = rhotc[cen_offset_0:x_a_dim[0] - cen_offset_0, cen_offset_1:x_a_dim[1] - cen_offset_1]
                gemo.write(bands[b]['rhot_ds'].replace('rhot', 'rhotc'), rhotc, ds_att = att)
                ## write also as rhot to L1RC file
                if setu['radcor_write_rhotc_separate_file']:
                    gemo_l1rc.write(bands[b]['rhot_ds'], rhotc, ds_att = att)

            ## write Ed
            if setu['output_ed']:
                ## baseline Ed with no spherical albedo effect
                if 'se_distance' not in gemo.gatts:
                    clon = np.nanmedian(gem.data('lon'))
                    clat = np.nanmedian(gem.data('lat'))
                    spos = ac.shared.sun_position(gemo.gatts['isodate'], clon, clat)
                    gemo.gatts['se_distance'] = spos['distance'] * 1.0
                Ed_base = (bands[b]['F0']) * gemo.gatts['se_distance']**2 * cos_sza * bands[b]['td_gas'] * T_d_tot
                att['Ed_base'] = Ed_base
                ## Ed with heterogeneous surface
                Ed =  Ed_base / (1 - rho_env_sph_est * rho_a_sph)
                if setu['radcor_crop_centre']: ## crop to centre area
                    Ed = Ed[cen_offset_0:x_a_dim[0] - cen_offset_0, cen_offset_1:x_a_dim[1] - cen_offset_1]
                gemo.write(bands[b]['rhot_ds'].replace('rhot', 'Ed'), Ed, ds_att = att)
                ## Ed with homogeneous surface
                #Ed = Ed_base / (1 - rho_s_homo * rho_a_sph)
                #if setu['radcor_crop_centre']: ## crop to centre area
                #    Ed = Ed[cen_offset_0:x_a_dim[0] - cen_offset_0, cen_offset_1:x_a_dim[1] - cen_offset_1]
                #gemo.write(bands[b]['rhot_ds'].replace('rhot', 'Edu'), Ed, ds_att = att)
                ## Ed with no surface - doesn't need to be written as array, but convenient for SNAP spectrum viewer
                #Ed[:] = Ed_base
                #gemo.write(bands[b]['rhot_ds'].replace('rhot', 'Ed0'), Ed, ds_att = att)

            ## write rhoe
            if setu['radcor_write_rhoe']:
                if setu['radcor_crop_centre']: ## crop to centre area
                    rho_env_psf = rho_env_psf[cen_offset_0:x_a_dim[0] - cen_offset_0, cen_offset_1:x_a_dim[1] - cen_offset_1]
                gemo.write(bands[b]['rhos_ds'].replace('rhos', 'rhoe'), rho_env_psf, ds_att = att)

            ## write rhos
            if setu['radcor_crop_centre']: ## crop to centre area
                rho_s_est = rho_s_est[cen_offset_0:x_a_dim[0] - cen_offset_0, cen_offset_1:x_a_dim[1] - cen_offset_1]
            gemo.write(bands[b]['rhos_ds'], rho_s_est, ds_att = att)

            ## write rho_s_homo
            if setu['radcor_write_rhosu']:
                if setu['radcor_crop_centre']: ## crop to centre area
                    rho_s_homo = rho_s_homo[cen_offset_0:x_a_dim[0] - cen_offset_0, cen_offset_1:x_a_dim[1] - cen_offset_1]
                gemo.write(bands[b]['rhos_ds'].replace('rhos_', 'rhosu_'), rho_s_homo, ds_att = att)

    ### end of inner correct_band function


    # --------------------------------------------------------------
    # Get scene metadata:
    #

    print('Running RAdCor on scene {}\n'.format(ncf))

    ## Open inputfile
    gem = ac.gem.gem(ncf)
    sensor = gem.gatts['sensor']

    ## check if L1R file is passed
    if gem.gatts['acolite_file_type'] != 'L1R':
        print('RAdCor processing not supported for acolite_file_type={}'.format(gem.gatts['acolite_file_type']))
        gem = None
        return

    ## combine default and user defined settings
    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}
    ## get sensor specific defaults
    setd = ac.acolite.settings.parse(sensor)
    ## set sensor default if user has not specified the setting
    for k in setd:
        if k not in ac.settings['user']: setu[k] = setd[k]
    ## end set sensor specific defaults
    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    ## Set up output directory
    output = setu['output']
    if output is None: output = os.path.dirname(ncf)

    ## File basename
    bn = os.path.basename(ncf)

    ## Set up oname and ofile
    oname = bn.replace('L1R', 'L2R').replace('.nc', '')
    if setu['region_name'] != '':
        if setu['region_name'] not in oname: oname+='_{}'.format(setu['region_name']) ## QV 2024-03-28
    ofile = '{}/{}.nc'.format(output, oname)

    if setu['rsr_version'] is not None:
        sensor_lut = '{}_{}'.format(sensor, setu['rsr_version'])
    else:
        sensor_lut = '{}'.format(sensor)

    ## if tsdsf_kernel_radius is not specified use radcor_kernel_radius
    if setu['tsdsf_kernel_radius'] is None: setu['tsdsf_kernel_radius'] = 1.0 * setu['radcor_kernel_radius'] ## QV 20240312

    ## validate radcor settings
    if not ac.adjacency.radcor.tools.validate_settings(setu):
        print('Error: Invalid settings for RAdCor... stopping.')
        gem = None
        return

    romix_par = 'romix'
    if setu['dsf_interface_reflectance']: romix_par = 'romix+rsky_t'
    add_rsky = romix_par == 'romix+rsky_t'
    ## additional imports
    #if setu['radcor_diagnostic_plots']:
    if True:
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable

    ## scene information
    ss = sensor.lower().split('_')
    x_a_dim = gem.ydim, gem.xdim # ac.adjacency.radcor.tools.nc_dim(ncf)
    print('Scene from {}.'.format(sensor))
    if sensor_lut != sensor: print('Using LUT for {}.'.format(sensor_lut))
    print('Scene dimensions: {}x{} (pixels)'.format(x_a_dim[0], x_a_dim[1]))

    ## START DEVELOPMENT BLOCK ##
    if setu['radcor_development']:
        dev_print = {}
        print('(DEV) Running RAdCor with:')
        for k in setu:
            if 'radcor' in k:
                print('    {} = {}'.format(k, setu[k]))
                dev_print[k] = setu[k]
        dev_print['file'] = ncf
        print('(DEV) {} = {}'.format('dsf_spectrum_option', setu['dsf_spectrum_option']))
    ## END DEVELOPMENT BLOCK ##

    ## Use scene_pixel_size attribute if available or get sensor defaults
    resolution = None
    if 'scene_pixel_size' in gem.gatts:
        resolution = gem.gatts['scene_pixel_size'][0]
        print('Got spatial resolution of {} m from global attributes.\n'.format(resolution))
    else:
        if sensor in ['S2A_MSI', 'S2B_MSI', 'S2C_MSI']:
            resolution = 10 ## can be 20 or 60 as well
            if 's2_target_res' in gem.gatts:
                resolution = gem.gatts['s2_target_res']
                print('Using resolution s2_target_res={} for MSI'.format(resolution))
        elif sensor in ['L8_OLI', 'L9_OLI']:
            resolution = 30 # 'PRISMA'
        elif sensor in ['EN1_MERIS', 'S3A_OLCI', 'S3B_OLCI']:
            if 's3_product_type' in gem.gatts:
                resolution = {'FR': 300, 'RR': 1200}[gem.gatts['s3_product_type']]
                print('Using resolution {} m based on s3_product_type={}'.format(resolution, gem.gatts['s3_product_type']))
            else:
                resolution = 300
                print('Using default resolution {} m'.format(resolution))
            print('Warning: Experimental RAdCor processing for {}'.format(sensor))
        elif 'PlanetScope' in sensor:
            resolution = 3
        elif sensor in ['PHR1A', 'PHR1B']:
            resolution = 2
            print('Warning: Experimental RAdCor processing for {}'.format(sensor))
        elif sensor in ['WorldView2', 'WorldView3']:
            resolution = 2
            print('Warning: Experimental RAdCor processing for {}'.format(sensor))
        elif sensor in ['WV_LG01', 'WV_LG01']:
            resolution = 1.36
            print('Warning: Experimental RAdCor processing for {}'.format(sensor))
        else:
            print('RAdCor processing not implemented for {}'.format(sensor))
            print('Not running RAdCor on scene {}'.format(ncf))
            gem = None
            return
        print('Will use {} sensor default spatial resolution: {}.\n'.format(sensor, resolution))

    if setu['dsf_interface_reflectance']:
        print('Warning: dsf_interface_reflectance=True, it is recommended to set it to False for RAdCor processing')

    ## Compute PSF radius in pixels
    psf_radius_pixels = int(setu['radcor_kernel_radius'] * 1000 / resolution) + 0.5

    ## Check PSF extent considering image extent
    if (x_a_dim[0] < psf_radius_pixels * 2) | (x_a_dim[1] < psf_radius_pixels * 2):
        print('Scene too small, it should cover at least 2 x PSF radius (2 x {} pixels)'.format(psf_radius_pixels))
        gem = None
        return

    ## aerosol models and names
    aer_nm = {'C': 'MOD1', 'M': 'MOD2'}
    aer_models = ["C", "M"]
#    aer_nm = {'C': 'MOD1', 'M': 'MOD2', 'U': 'MOD3'} # ALEX 2024-08-27
#    aer_models = ["C", "M", "U"] # ALEX 20240827
    if setu['radcor_force_model'] is not None: ## do only forced model
        aer_models = [setu['radcor_force_model'].upper()]
        if aer_models[0] not in aer_nm:
            print('Model radcor_force_model={} not configured'.format(aer_models[0]))
            return

    ## Get scene average geometry
    sza, vza, raa = gem.gatts['sza'], gem.gatts['vza'], gem.gatts['raa']
    if vza > setu['radcor_max_vza']:
        print('The current implementation of RAdCor assumes a circular PSF and is not suited for higher viewing zenith angles.')
        print('Scene average viewing zenith angle {:.1f} > radcor_max_vza={:.1f}'.format(np.nanmean(vza), setu['radcor_max_vza']))
        return

    cos_sza = np.cos(np.radians(sza))
    cos_vza = np.cos(np.radians(vza))

    ## Load RSR dict
    if sensor in ac.config['hyper_sensors']:
        rsr = ac.shared.rsr_hyper(gem.gatts['band_waves'], gem.gatts['band_widths'], step=0.1)
        rsrd = ac.shared.rsr_dict(rsrd = {sensor_lut : {'rsr' : rsr}})
    else:
        rsrd = ac.shared.rsr_dict(sensor = sensor_lut)

    #
    # End Get scene metadata
    # --------------------------------------------------------------

    # --------------------------------------------------------------
    # Get ancillary data:
    #

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

    ## get gas transmittance
    tg_dict = ac.ac.gas_transmittance(sza, vza, uoz = uoz, uwv = uwv, rsr = rsrd[sensor_lut]['rsr'])

    ## get F0 and downward gas transmittance to compute Ed
    if setu['output_ed']:
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data']), rsrd[sensor_lut]['rsr'])
        ## downward gas transmittance
        tdgas = ac.ac.gaslut_interp(sza, vza, pressure=pressure,
                                    sensor=None, lutconfig='202402F', pars = ['dtdica','dtoxyg','dtniox','dtmeth'])
        tdg = np.prod([tdgas[k] for k in tdgas if k != 'wave'], axis=0)
        tdg_b = ac.shared.rsr_convolute_dict(np.asarray(tdgas['wave']), tdg, rsrd[sensor_lut]['rsr'])

    ## START DEVELOPMENT BLOCK ##
    if setu['radcor_development']:
        dev_print['sensor'] = sensor
        dev_print['resolution'] = resolution
        dev_print['dimensions'] = x_a_dim
        print('(DEV) Sun zenith angle: {} (degrees)'.format(sza))
        print('(DEV) View zenith angle: {} (degrees)'.format(vza))
        print('(DEV) Relative azimuth angle: {} (degrees)'.format(raa))
        dev_print['sza'] = sza
        dev_print['vza'] = vza
        dev_print['raa'] = raa
        dev_print['pressure'] = pressure
        dev_print['tt_gas_'] = tg_dict['tt_gas']
        print('(DEV) Ozone column: {:.3f}'.format(uoz))
        print('(DEV) Water vapour: {:.3f}'.format(uwv))
        print('(DEV) Atmospheric pressure at the surface: {:.2f}'.format(pressure))
        print('\n(DEV) Gas transmittance:')
        for bi, b in enumerate(rsrd[sensor_lut]['rsr_bands']):
            print('    Band {}: {:.5f}'.format(b, tg_dict['tt_gas'][b]))
    ## END DEVELOPMENT BLOCK ##

    #
    # End Get ancillary data
    # --------------------------------------------------------------


    # --------------------------------------------------------------
    # Create bands dictionary and select bands:
    #

    ## Bands dictionary
    bands = {}
    bint = 0
    for bi, b in enumerate(rsrd[sensor_lut]['rsr_bands']):
        if b not in bands:
            rhot_ds = 'rhot_{}'.format(rsrd[sensor_lut]['wave_name'][b])
            if rhot_ds not in gem.datasets:
                print('{} dataset for band {} not in inputfile.'.format(rhot_ds, b))
                continue

            bands[b] = {k:rsrd[sensor_lut][k][b] for k in ['wave_mu', 'wave_nm', 'wave_name'] if b in rsrd[sensor_lut][k]}
            bands[b]['rhot_ds'] = 'rhot_{}'.format(bands[b]['wave_name'])
            bands[b]['rhos_ds'] = 'rhos_{}'.format(bands[b]['wave_name'])
            for k in tg_dict:
                if k not in ['wave']: bands[b][k] = tg_dict[k][b]
            bands[b]['wavelength'] = bands[b]['wave_nm']
            if setu['output_ed']:
                bands[b]['F0'] = f0_b[b]
                bands[b]['td_gas'] = tdg_b[b]
            ## track APSFS simulation band naming
            bint += 1
            ## subset SD8 bands for SD5
            if sensor == 'PlanetScope_SD5':
                aliases = {1:2, 2:4, 3:6, 4:7, 5:8}
                bands[b]['fit_name'] = aliases[bint]
            # elif sensor == 'L5_TM':
            #     aliases = {1:1, 2:2, 3:3, 4:4, 5:5, 7:6}
            #     bands[b]['fit_name'] = aliases[bint]
            # elif sensor == 'L7_ETM':
            #     aliases = {1:1, 2:2, 3:3, 4:4, 5:5, 7:6, 8:7}
            #     bands[b]['fit_name'] = aliases[bint]
            else:
                bands[b]['fit_name'] = bint ## band name (number) in fits files
            ## track use of band in RAdCor/TSDSF
            bands[b]['radcor_use_band'] = True
            bands[b]['tsdsf_use_band'] = True
            if bands[b]['tt_gas'] < setu['min_tgas_rho']: bands[b]['radcor_use_band'] = False
            if bands[b]['tt_gas'] < setu['min_tgas_aot']: bands[b]['radcor_use_band'] = False
            if (setu['radcor_skip_pan']) & ('OLI' in sensor) & (b == '8'): bands[b]['radcor_use_band'] = False
            if (bands[b]['wavelength'] < setu['tsdsf_wave_range'][0]): bands[b]['tsdsf_use_band'] = False
            if (bands[b]['wavelength'] > setu['tsdsf_wave_range'][1]): bands[b]['tsdsf_use_band'] = False

    ## List of bands to be used
    bands_ = [b for b in bands if (bands[b]['radcor_use_band']) & (bands[b]['tsdsf_use_band'])]
    if len(bands_) == 0:
        print('No band data found in inputfile.')
        gem = None
        return

    ## START DEVELOPMENT BLOCK ##
    if setu['radcor_development']:
        print('(DEV) Selected bands: {}'.format(bands_))
        dev_print['selected_bands'] = bands_
        dev_print['lambda'] = [bands[b]['wavelength'] for b in bands_]
    ## END DEVELOPMENT BLOCK ##

    #
    # End Create bands dictionary and select bands
    # --------------------------------------------------------------

    # --------------------------------------------------------------
    ## Test if settings for aot optimisation to field measurement are correct
    #
    if setu['radcor_aot_estimate'] == 'optimise':
        print()

        ## if optimise_target_rhos_file is given, read the data, and resample to the sensor RSR
        if (setu['optimise_target_rhos_file'] is not None):
            ## read file for optimisation
            ret = ac.ac.optimise_read(setu)
            if ret is None:
                gem = None
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
            gem = None
            return

        ## check if number of target rhos corresponds to the number of considered bands
        if len(setu['optimise_target_rhos']) != len(bands_):
            print('The number of items in optimise_target_rhos ({}) does not match the number of considered bands ({}).'.format(len(setu['optimise_target_rhos']), len(bands_)))
            print('Provide optimise_target_rhos for each considered band: {}'.format(bands_))
            print('Set missing bands (e.g. SWIR) to NaN to be ignored, or to 0 to take them into account in the fit')
            gem = None
            return

        ## don't provide all NaNs!
        if not any(np.isfinite(np.asarray(setu['optimise_target_rhos'], dtype=float))):
            print('Zero finite items in optimise_target_rhos: {}'.format(', '.join([str(v) for v in setu['optimise_target_rhos']])))
            gem = None
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
            gem = None
            return
    #
    ## End test aot optimisation to field measurement settings
    # --------------------------------------------------------------


    # --------------------------------------------------------------
    # Import LUTs:
    #

    lutdw = ac.aerlut.import_luts(add_rsky = add_rsky, par = romix_par, sensor = sensor_lut,
        lut_par = ['utott', 'dtott', 'astot', 'romix', 'tray', 'taer', 'ttot', 'dutott', 'utotr', 'utota', 'asray', 'asaer'], ## AC 2024-09-09 ## Added asray and asaer for SAF calculation
#        base_luts = ['ACOLITE-LUT-202110-MOD1', 'ACOLITE-LUT-202110-MOD2', 'ACOLITE-LUT-202110-MOD3'], # ALEX 2024-08-27
        #return_lut_array = True)
        return_lut_array = False)

    luts = list(lutdw.keys())

    #
    # End Import LUTs
    # --------------------------------------------------------------


    # --------------------------------------------------------------
    # APSFs and SAFs data:
    #

    ## import APSF and SAF data from separate function
    coefs_psf_ray, coefs_psf_aer, coefs_saf_ray, coefs_saf_aer = ac.adjacency.radcor.tools.import_fits(sensor, bands, aer_models = aer_models)

    ## START DEVELOPMENT BLOCK ##
    if setu['radcor_development']:
        for am in coefs_psf_aer:
            print('(DEV) PSF model {}:'.format(am))
            print('    Pressure dependence: {}'.format(coefs_psf_aer[am]['1']['press_dep']))
            print('    Pressure Range: {}'.format(coefs_psf_aer[am]['1']['press_rng']))
            print('    Coefficients:')
            for k in coefs_psf_aer[am]['1']['coefficients']:
                print('        {}: {:.3f}'.format(k, coefs_psf_aer[am]['1']['coefficients'][k]))
            print('\n')

        am = 'fa_generic_R'
        print('(DEV) PSF model R:')
        print('    Pressure dependence: {}'.format(coefs_psf_ray[am]['press_dep']))
        print('    Pressure Range: {}'.format(coefs_psf_ray[am]['press_rng']))
        print('    Coefficients:')
        for k in coefs_psf_ray[am]['coefficients']:
            print('        {}: {:.3f}'.format(k, coefs_psf_ray[am]['coefficients'][k]))
        print('\n\n')

        for am in coefs_saf_aer:
            print('(DEV) SAF model {}:'.format(am))
            print('    Pressure dependence: {}'.format(coefs_saf_aer[am]['1']['press_dep']))
            print('    Pressure Range: {}'.format(coefs_saf_aer[am]['1']['press_rng']))
            print('    Coefficients:')
            for k in coefs_saf_aer[am]['1']['coefficients']:
                print('        {}: {:.3f}'.format(k, coefs_saf_aer[am]['1']['coefficients'][k]))
            print('\n')

        am = 'fa_generic_R'
        print('(DEV) SAF model R:')
        print('    Pressure dependence: {}'.format(coefs_saf_ray[am]['press_dep']))
        print('    Pressure Range: {}'.format(coefs_saf_ray[am]['press_rng']))
        print('    Coefficients:')
        for k in coefs_saf_ray[am]['coefficients']:
            print('        {}: {:.3f}'.format(k, coefs_saf_ray[am]['coefficients'][k]))
        print('\n\n')
    ## END DEVELOPMENT BLOCK ##

    #
    # End APSFs and SAFs data
    # --------------------------------------------------------------

    # --------------------------------------------------------------
    # Run AOT estimation
    #

    ## Track rho_dark and rho_toa_av per band
    nbands_    = len(bands_)
    rho_dark   = np.zeros(nbands_)
    rho_toa_av = np.zeros(nbands_)

    ## Track rho_a and AOT estimation per aerosol model and per band
    rho_a_est  = np.zeros((3, nbands_)) #ALEX 2024-08-27
    tau550_est = np.zeros((3, nbands_)) #ALEX 2024-08-27
    idmin      = np.zeros((3, nbands_, 2)) #ALEX 2024-08-27

    ## user fixed
    if (setu['radcor_force_model'] is not None) & (setu['radcor_force_aot'] is not None):
        print('\nUser supplied model and aot: radcor_force_model={} radcor_force_aot={}'.format(setu['radcor_force_model'],setu['radcor_force_aot']))
        best_mod = setu['radcor_force_model']
        best_aot = setu['radcor_force_aot']
        ## set these to nan if user forces model and aot
        best_fit = np.nan
        best_idx = np.nan
        best_band = np.nan
    else:
        print('\nRunning AOT estimation using radcor_aot_estimate={}'.format(setu['radcor_aot_estimate']))

        if setu['radcor_aot_estimate'] == 'dsf':
            # --------------------------------------------------------------
            # Run DSF:
            #

            ## offsets to extract centre of the image - for comparison with TS-DSF
            cen_offset_0 = 0
            cen_offset_1 = 0
            ## change offsets if radcor_aot_estimate_centre_extent is given
            if setu['radcor_aot_estimate_centre_extent'] is not None:
                cen_offset_npix = np.round(float(setu['radcor_aot_estimate_centre_extent'] * 1000) / resolution).astype(int)
                cen_offset_0 = np.round((x_a_dim[0] - cen_offset_npix) / 2).astype(int)
                cen_offset_1 = np.round((x_a_dim[1] - cen_offset_npix) / 2).astype(int)
                print('Limiting scene centre to radcor_aot_estimate_centre_extent={:.1f} km ({} pixels)'.format(setu['radcor_aot_estimate_centre_extent'], cen_offset_npix))
                print('Offsets to extract image centre: {}x{} pixels'.format(cen_offset_0, cen_offset_1))
                print('Extracted centre size: {}x{} pixels'.format(x_a_dim[0] - 2*cen_offset_0, x_a_dim[1] - 2 * cen_offset_1))
            ## end offsets to extract centre of the image

            ## Run through bands to estimate AOT
            for bi, b in enumerate(bands_):
                ## Load rhot data
                print('\nLoading TOA data for band {}: {}'.format(b, bands[b]['rhot_ds']))
                rho_toa, att = gem.data(bands[b]['rhot_ds'], attributes = True)

                ## crop to centre (if offsets != 0)
                if setu['radcor_aot_estimate_centre_extent'] is not None:
                    rho_toa = rho_toa[cen_offset_0:x_a_dim[0] - cen_offset_0, cen_offset_1:x_a_dim[1] - cen_offset_1] ## QV 2024-03-26
                    print(rho_toa.shape)

                ## Gas transmittance correction
                rho_toa /= bands[b]['tt_gas']

                ## Scene stats
                rho_toa_av[bi] = np.nanmean(rho_toa)

                ## Get rho_dark in whole scene
                if setu['dsf_spectrum_option'] == 'darkest':
                    rho_dark[bi] = np.array((np.nanpercentile(rho_toa, 0)))
                if setu['dsf_spectrum_option'] == 'percentile':
                    rho_dark[bi] = np.array((np.nanpercentile(rho_toa, setu['dsf_percentile'])))
                if setu['dsf_spectrum_option'] == 'intercept':
                    rho_dark[bi] = ac.shared.intercept(rho_toa, setu['dsf_intercept_pixels'])

                print(setu['dsf_spectrum_option'], setu['dsf_intercept_pixels'])
                ## Estimate AOT(550) for each model
                print('    Estimating aot550 for models {}'.format(', '.join(aer_models)))
                for ai, am in enumerate(aer_models):
                    lut = [lut for lut in luts if aer_nm[am] in lut][0]
                    ## path reflectance according to tau_a in LUT
                    if add_rsky:
                        romix_ = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd'][romix_par], raa, vza, sza, wind, lutdw[lut]['meta']['tau']))
                    else:
                        romix_ = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['romix'], raa, vza, sza, lutdw[lut]['meta']['tau']))
                    ## find tau according to observed rho_dark
                    tau550_est[ai, bi] = np.interp(rho_dark[bi], romix_, lutdw[lut]['meta']['tau'])
                    rho_a_est[ai, bi] = 1.0 * rho_dark[bi] ## set here to observed rho_dark

                    print('    band {} MOD {} rho_dark {:.4f} tau_a 550 nm {:.4f}'.format(b, am, rho_dark[bi], tau550_est[ai, bi]))
            #
            # End DSF
            # --------------------------------------------------------------

        elif setu['radcor_aot_estimate'] == 'tsdsf':
            # --------------------------------------------------------------
            # Run TS-DSF:
            #

            ## x vector for TS-DSF PSF:
            xvec = np.arange(resolution / 2000.0, setu['tsdsf_kernel_radius'], step = resolution / 1000.0)
            xvec = np.hstack((-1.0 * np.flip(xvec), xvec))[1:]
            id_psf = int((len(xvec) - 1) / 2)
            x_f_dim = x_a_dim * 1
            x_o_dim = x_a_dim * 1
            if not setu['radcor_edge_extend']:
                x_o_dim = int(x_a_dim[0] - id_psf * 2), int(x_a_dim[1] - id_psf * 2)
            else:
                x_f_dim = int(x_a_dim[0] + id_psf * 2), int(x_a_dim[1] + id_psf * 2)

            ## offsets to extract centre of the image
            cen_offset_0 = int(id_psf)
            cen_offset_1 = int(id_psf)
            ## change offsets if radcor_aot_estimate_centre_extent is given
            if setu['radcor_aot_estimate_centre_extent'] is not None:
                cen_offset_npix = np.round(float(setu['radcor_aot_estimate_centre_extent'] * 1000) / resolution).astype(int)
                cen_offset_0 = np.round((x_a_dim[0] - cen_offset_npix) / 2).astype(int)
                cen_offset_1 = np.round((x_a_dim[1] - cen_offset_npix) / 2).astype(int)
                print('Limiting TS-DSF scene centre to radcor_aot_estimate_centre_extent={:.1f} km ({} pixels)'.format(setu['radcor_aot_estimate_centre_extent'], cen_offset_npix))
            else:
                print('Limiting TS-DSF scene to PSF coverage with radius of {:.1f} km ({} pixels)'.format(setu['radcor_kernel_radius'], psf_radius_pixels))
            print('Offsets to extract image centre: {}x{} pixels'.format(cen_offset_0, cen_offset_1))
            print('Extracted centre size: {}x{} pixels'.format(x_a_dim[0] - 2*cen_offset_0, x_a_dim[1] - 2 * cen_offset_1))
            ## end offsets to extract centre of the image

            ## Create uniform PSF and OTF
            if (setu['tsdsf_kernel_complete_method'] == "neighbourhood"):
                psf_uni = np.zeros((len(xvec), len(xvec))) + 1.0 / len(xvec)**2
                #otf_uni = np.zeros(x_f_dim) ## AC 2024-09-09 ## No longer necessary as psf_uni is summed to psf_mix and saf_mix
                #otf_uni[0:psf_uni.shape[0], 0:psf_uni.shape[1]] = psf_uni ## AC 2024-09-09 ##
                #otf_uni = np.fft.fftn(otf_uni) ## AC 2024-09-09 ##

            ## START DEVELOPMENT BLOCK ##
            if setu['radcor_development']:
                print('Length of xvec = {}'.format(len(xvec)))
                if (setu['tsdsf_kernel_complete_method'] == "neighbourhood"):    # AC 20240306
                    print('Value of PSF uni = {}'.format(psf_uni[0][0]))
                    print('Integral of psf_uni = {}'.format(np.sum(psf_uni)))

                psf_tsdsf_mix_C_cov = np.zeros(nbands_)
                psf_tsdsf_mix_M_cov = np.zeros(nbands_)

                psf_ofile = '{}/{}'.format(output, bn.replace('L1R', '').replace('.nc', '_psf.nc'))
                gempsf = ac.gem.gem(psf_ofile, new = True)
                gempsf.gatts = {k: gem.gatts[k] for k in gem.gatts}
                gempsf.gatts['auto_grouping'] = 'rhot:rhorc:rhos:rhow:Rrs:rhoe:rho_nbh:psf_uni:psf_tsdsf_mix_C:psf_tsdsf_mix_M:psf_mix:psf_scl'
                print(psf_ofile)
                if (setu['tsdsf_kernel_complete_method'] == "neighbourhood"):    # AC 20240306
                    a_ = {'units': '1', '_FillValue': 1e+30, 'long_name': 'Discrete normalized uniform PSF', 'parameter': 'psf_uni', 'standard_name': 'normalized_uniform_psf'}
                    gempsf.write('psf_uni', psf_uni, ds_att = a_)
            ## END DEVELOPMENT BLOCK ##

            ## Run through bands to estimate AOT
            for bi, b in enumerate(bands_):
                ## Load rhot data
                print('\nLoading TOA data for band {}: {}'.format(b, bands[b]['rhot_ds']))
                rho_toa, att = gem.data(bands[b]['rhot_ds'], attributes = True)

                ## Get Rayleigh optical depth
                if add_rsky:
                    tau_ray = lutdw[luts[0]]['rgi'][b]((pressure, lutdw[luts[0]]['ipd']['tray'], raa, vza, sza, wind, 0.001))
                else:
                    tau_ray = lutdw[luts[0]]['rgi'][b]((pressure, lutdw[luts[0]]['ipd']['tray'], raa, vza, sza, 0.001))
                print('    Rayleigh optical depth: {:.3f} (unitless)'.format(tau_ray))

                ## Gas transmittance correction
                rho_toa /= bands[b]['tt_gas']

                ## Scene stats
                rho_toa_av[bi] = np.nanmean(rho_toa)

                ## Array of rho_toa
                x_a = 1.0 * rho_toa

                ## Fill NaNs
                x_a[np.isnan(x_a)] = rho_toa_av[bi]

                ## Extend rhot array if requested to keep edges
                if setu['radcor_edge_extend']:
                    x_a = ac.adjacency.radcor.extend.extend(x_a, id_psf, x_a_dim, x_f_dim, method = setu['radcor_edge_extend_method'])

                ## Do FFT of extended TOA data
                print('    Computing FFT of TOA array')
                x_f = np.fft.fftn(x_a)

                ## Estimate neighbourhood
                #if (not setu['tsdsf_kernel_rescale']) & (setu['tsdsf_kernel_complete_method'] == "neighbourhood"): ## AC 2024-09-09 ## No longer necessary as psf_uni is summed to psf_mix
                #    rho_nbh_toa = ac.adjacency.radcor.tools.conv(x_f, otf_uni, id_psf, x_f_dim, keep_edges = setu['radcor_edge_extend']) ## AC 2024-09-09 ##

                ## START DEVELOPMENT BLOCK ##
                if setu['radcor_development']:
                    if 'conv_ofile' not in locals():
                        conv_ofile = '{}/{}'.format(output, bn.replace('L1R', '').replace('.nc', '_conv.nc'))
                        gemconv = ac.gem.gem(conv_ofile, new = True)
                        gemconv.gatts = {k: gempsf.gatts[k] for k in gempsf.gatts}
                    #if (not setu['tsdsf_kernel_rescale']) & (setu['tsdsf_kernel_complete_method'] == "neighbourhood"):    # AC 20240306 ## AC 2024-09-09 ## Removed because this average is no longer used.
                    #    a_ = {'units': '1', '_FillValue': 1e+30, 'long_name': 'TOA uniform neighbourhood average',
                    #        'parameter': 'rho_nbh_' + bands[b]['wave_name'], 'standard_name': 'uniform_neighbourhood_average',
                    #        'wavelength': bands[b]['wave_name']} ## AC 2024-09-09 ##
                    #    gemconv.write('rho_nbh_' + bands[b]['wave_name'], rho_nbh_toa, ds_att = a_) ## AC 2024-09-09 ##
                ## END DEVELOPMENT BLOCK ##

                ## Crop to centre of scene fully covered by PSF or to requested centre
                #rho_toa_cen = rho_toa[id_psf:x_a_dim[0] - id_psf, id_psf:x_a_dim[1] - id_psf]
                rho_toa_cen = rho_toa[cen_offset_0:x_a_dim[0] - cen_offset_0, cen_offset_1:x_a_dim[1] - cen_offset_1] ## QV 2024-03-26

                ## Get rho_dark in this central portion:
                if setu['dsf_spectrum_option'] == 'darkest':
                    rho_dark[bi] = np.array((np.nanpercentile(rho_toa_cen, 0)))
                if setu['dsf_spectrum_option'] == 'percentile':
                    rho_dark[bi] = np.array((np.nanpercentile(rho_toa_cen, setu['dsf_percentile'])))
                if setu['dsf_spectrum_option'] == 'intercept':
                    rho_dark[bi] = ac.shared.intercept(rho_toa_cen, setu['dsf_intercept_pixels'])

                print('    Scene shape {}'.format(rho_toa.shape))
                print('    Scene centre shape {}'.format(rho_toa_cen.shape))

                ## START DEVELOPMENT BLOCK ##
                if setu['radcor_development']:
                    print('\n    (DEV) Average rho_toa: {:.3f}'.format(rho_toa_av[bi]))
                    print('    (DEV) Range of rho_toa: {:.3f}'.format(rho_dark[bi], np.nanmax(rho_toa)))
                    print('    (DEV) Range of x_a: {:.3f} - {:.3f}'.format(np.nanmin(x_a), np.nanmax(x_a)))
                    print('    (DEV) Range of real part of FFT: {:.3f} - {:.3f}'.format(np.nanmin(np.real(x_f)), np.nanmax(np.real(x_f))))
                    #if (not setu['radcor_kernel_rescale']) & (setu['radcor_kernel_complete_method'] == "neighbourhood"): # AC 20240306
                    #if (not setu['tsdsf_kernel_rescale']) & (setu['tsdsf_kernel_complete_method'] == "neighbourhood"):    # AC 20240306 ## AC 2024-09-09 ## psf_uni no longer separately used
                    #    print('    (DEV) Uniform PSF value: {:.7f} (unitless)\n'.format(psf_uni[0][0]))
                ## END DEVELOPMENT BLOCK ##

                ## Estimate AOT(550) for each model
                print('    Estimating aot550 for models {}'.format(', '.join(aer_models)))
                for ai, am in enumerate(aer_models):
                    lut = [lut for lut in luts if aer_nm[am] in lut][0]

                    ## Aerosol coefficients for sensor bands from previously imported data
                    aer_psf = coefs_psf_aer[am][b]

                    if add_rsky:
                        T_u_ray = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['utotr'], raa, vza, sza, wind, setu['tsdsf_initial_aot']))
                        T_u_aer = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['utota'], raa, vza, sza, wind, setu['tsdsf_initial_aot']))
                        tau_aer = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['taer'], raa, vza, sza, wind, setu['tsdsf_initial_aot']))
                    else:
                        T_u_ray = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['utotr'], raa, vza, sza, setu['tsdsf_initial_aot']))
                        T_u_aer = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['utota'], raa, vza, sza, setu['tsdsf_initial_aot']))
                        tau_aer = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['taer'], raa, vza, sza, setu['tsdsf_initial_aot']))
                    T_u_dif_ray = T_u_ray - np.exp(-tau_ray / cos_vza)
                    T_u_dif_aer = T_u_aer - np.exp(-tau_aer / cos_vza)

                    coef_psf_aer = [aer_psf['coefficients'][c] for c in aer_psf['coefficients']]
                    coef_psf_ray = [coefs_psf_ray['fa_generic_R']['coefficients'][c] for c in coefs_psf_ray['fa_generic_R']['coefficients']]

                    ## Set up PSF kernel
                    psf_mix, ext = ac.adjacency.radcor.tools.w_kernel(coef_psf_aer, coef_psf_ray, T_u_dif_ray, T_u_dif_aer,
                        ex = setu['tsdsf_kernel_radius'], res = resolution / 1000, pressure = pressure)
                    psf_cv_mix = np.sum(psf_mix)
                    print('    Model {} T_u_dif_ray {:.4f}'.format(am, T_u_dif_ray))
                    print('    Model {} T_u_dif_aer {:.4f}'.format(am, T_u_dif_aer))
                    print('    Model {} PSF cover {:.4f}'.format(am, psf_cv_mix))
                    if (setu['radcor_kernel_complete_method'] == "renormalise"):
                        psf_mix /= psf_cv_mix
                        psf_cv_mix = 1.0
                    elif (setu['radcor_kernel_complete_method'] == "neighbourhood"): ## AC 2024-09-09 ## psf_uni now added directly to the psf_mix
                        psf_mix += (1.0 - psf_cv_mix) * psf_uni ## AC 2024-09-09 ##
                        psf_cv_mix = 1.0 ## AC 2024-09-09 ##

                    ## test for setting TS-DSF PSF centre pixel to zero
                    if setu['tsdsf_psf_centre_zero']:
                        print('    Setting psf_mix centre pixel to zero')
                        psf_mix[id_psf, id_psf] = 0.0

                    ## START DEVELOPMENT BLOCK ##
                    if setu['radcor_development']:
                        print('Mix PSF shape {}'.format(psf_mix.shape))
                        #if (not setu['tsdsf_kernel_rescale']) & (setu['tsdsf_kernel_complete_method'] == "neighbourhood"):    # AC 20240306 ## AC 2024-09-09 ## psf_uni no longer separately used
                        #    print('Uni PSF shape {}'.format(psf_uni.shape)) ## AC 2024-09-09 ##

                        a_ = {'units': '1', '_FillValue': 1e+30, 'long_name': 'Discrete normalized mixed diffuse PSF',
                              'parameter': 'psf_tsdsf_mix_' + am + '_' + bands[b]['wave_name'], 'standard_name': 'normalized_mixed_diffuse_psf',
                              'wavelength': bands[b]['wave_name'], 'coverage': psf_cv_mix}
                        gempsf.write('psf_tsdsf_mix_' + am + '_' + bands[b]['wave_name'], psf_mix, ds_att = a_)

                        if am == 'C':
                            psf_tsdsf_mix_C_cov[bi] = psf_cv_mix
                        if am == 'M':
                            psf_tsdsf_mix_M_cov[bi] = psf_cv_mix
                    ## END DEVELOPMENT BLOCK ##

                    ## Compute OTF of mix PSF
                    otf_mix = np.zeros(x_f_dim)
                    otf_mix[0:psf_mix.shape[0], 0:psf_mix.shape[1]] = 1.0 * psf_mix
                    otf_mix = np.fft.fftn(otf_mix)

                    ## Estimate rho_env at TOA
                    rho_env_toa_est = ac.adjacency.radcor.tools.conv(x_f, otf_mix, id_psf, x_f_dim, keep_edges = setu['radcor_edge_extend'])

                    if (setu['tsdsf_kernel_complete_method'] == "average"):
                        rho_env_toa_est += (1 - psf_cv_mix) * rho_toa_av[bi]
                    #if (setu['tsdsf_kernel_complete_method'] == "neighbourhood"): ## AC 2024-09-09 ## psf_uni was added to psf_mix
                    #    rho_env_toa_est += (1 - psf_cv_mix) * rho_nbh_toa

                    ## START DEVELOPMENT BLOCK ##
                    if setu['radcor_development']:
                        a_ = {'units': '1', '_FillValue': 1e+30, 'long_name': 'Hypothetical quantity estimate',
                              'parameter': 'rho_hyp_' + am + '_' + bands[b]['wave_name'], 'standard_name': 'hypothetical_quantity_estimate',
                              'wavelength': bands[b]['wave_name']}
                        gemconv.write('rho_hyp_' + am + '_' + bands[b]['wave_name'], rho_env_toa_est, ds_att = a_)
                    ## END DEVELOPMENT BLOCK ##

                    ## Crop to centre covered by PSF or to requested centre
                    #rho_env_toa_est_cen = rho_env_toa_est[id_psf:x_a_dim[0] - id_psf, id_psf:x_a_dim[1] - id_psf]
                    rho_env_toa_est_cen = rho_env_toa_est[cen_offset_0:x_a_dim[0] - cen_offset_0, cen_offset_1:x_a_dim[1] - cen_offset_1] ## QV 2024-03-26

                    ## Estimate AOT(550)
                    if (setu['radcor_force_aot'] is not None):
                        tau550_est[ai, bi] = setu['radcor_force_aot'] * 1.0
                        if add_rsky:
                            rho_a_est[ai, bi]  = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd'][romix_par],raa, vza, sza, wind, tau550_est[ai, bi]))
                        else:
                            rho_a_est[ai, bi]  = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['romix'],raa, vza, sza, tau550_est[ai, bi]))
                    else:
                        ## Ratio squared of rho_toa to estimated rho_env_toa
                        bratio = rho_toa_cen * rho_toa_cen / rho_env_toa_est_cen

                        ## LUT interpolated to observation geometry
                        if add_rsky:
                            lut_rho_a   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd'][romix_par], raa, vza, sza, wind, lutdw[lut]['meta']['tau']))
                            lut_T_d_tot = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['dtott'], raa, vza, sza, wind, lutdw[lut]['meta']['tau']))
                            lut_T_u_tot = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['utott'], raa, vza, sza, wind, lutdw[lut]['meta']['tau']))
                            lut_tau_tot = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['ttot'],  raa, vza, sza, wind, lutdw[lut]['meta']['tau']))
                        else:
                            lut_rho_a   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['romix'], raa, vza, sza, lutdw[lut]['meta']['tau']))
                            lut_T_d_tot = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['dtott'], raa, vza, sza, lutdw[lut]['meta']['tau']))
                            lut_T_u_tot = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['utott'], raa, vza, sza, lutdw[lut]['meta']['tau']))
                            lut_tau_tot = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['ttot'],  raa, vza, sza, lutdw[lut]['meta']['tau']))

                        lut_T_u_dir   = np.exp(-lut_tau_tot / cos_vza)
                        lut_T_u_dif   = lut_T_u_tot - lut_T_u_dir
                        lut_T_u_dif_r = lut_T_u_dif / lut_T_u_tot

                        rho_a_opt_mn = lut_rho_a[0]
                        rho_a_opt_mx = rho_dark[bi]

                        # This is a quick fix to avoid using this band if the minimum rho_a in the LUT is greater than rho_toa...
                        if rho_a_opt_mn > rho_a_opt_mx:
                            print('    Warning: minimum rho_a in LUT {:.5f} is greater than rho_toa {:.5f}'.format(rho_a_opt_mn, rho_a_opt_mx))
                            tau = 5.0
                            opt = np.interp(tau, lutdw[lut]['meta']['tau'], lut_rho_a)
                        else:
                            # Avoid using bad pixels for finding the darkest pixel at the surface
                            bratio[rho_toa_cen < rho_dark[bi]] = np.nan
                            id = np.where(np.isfinite(bratio))
                            if setu['tsdsf_bratio_option'] == 'darkest':
                                minbr = np.array((np.nanpercentile(bratio[id], 0, method = 'closest_observation')))
                            if setu['tsdsf_bratio_option'] == 'percentile':
                                minbr = np.array((np.nanpercentile(bratio[id], setu['tsdsf_bratio_percentile'], method = 'closest_observation')))
                            minidx = np.where(bratio == minbr)
                            idmin[ai,bi,:] = minidx[0][0], minidx[1][0]

                            ## Find AOT(550) by optimization
                            start = np.nanmean((rho_a_opt_mn, rho_a_opt_mx))
                            a_min = rho_toa_cen[minidx[0][0], minidx[1][0]]
                            b_min = rho_env_toa_est_cen[minidx[0][0], minidx[1][0]]

                            opt = scipy.optimize.minimize_scalar(ac.adjacency.radcor.tools.ts_dsf_optfun, args = (a_min, b_min, lut_rho_a, lut_T_u_dif_r), bounds = (rho_a_opt_mn, rho_a_opt_mx), method = 'bounded').x
                            tau = np.interp(opt, lut_rho_a, lutdw[lut]['meta']['tau'])

                        rho_a_est[ai,bi] = opt
                        tau550_est[ai,bi] = tau

                        ## START DEVELOPMENT BLOCK ##
                        if setu['radcor_diagnostic_plots']:
                            ## plot rho_toa
                            fig, ax = plt.subplots()
                            ax.set_title('rho_toa B{} MOD{}'.format(b, am));
                            im = ax.imshow(rho_toa)
                            divider = make_axes_locatable(ax)
                            cax = divider.append_axes("right", size = "5%", pad = 0.05)
                            plt.colorbar(im, cax = cax)
                            fname = 'B{}_MOD{}_{}'.format(b, am, 'rho_toa')
                            plt.savefig('{}/{}_{}.png'.format(output, oname, fname), dpi = 300, bbox_inches = 'tight', facecolor = 'white')
                            plt.close()

                            #if (not setu['tsdsf_kernel_rescale']) & (setu['tsdsf_kernel_complete_method'] == "neighbourhood"): ## AC 2024-09-09 ## psf_uni no longer separately used
                            #    ## plot rho_nbh_toa
                            #    fig, ax = plt.subplots()
                            #    ax.set_title('rho_nbh_toa B{} MOD{}'.format(b, am));
                            #    im = ax.imshow(rho_nbh_toa)
                            #    divider = make_axes_locatable(ax)
                            #    cax = divider.append_axes("right", size = "5%", pad = 0.05)
                            #    plt.colorbar(im, cax = cax)
                            #    fname = 'B{}_MOD{}_{}'.format(b, am, 'rho_nbh_toa')
                            #    plt.savefig('{}/{}_{}.png'.format(output, oname, fname), dpi = 300, bbox_inches = 'tight', facecolor = 'white')
                            #    plt.close()

                            ## plot rho_env_toa_est
                            fig, ax = plt.subplots()
                            ax.set_title('rho_env_toa_est B{} MOD{}'.format(b, am));
                            im = ax.imshow(rho_env_toa_est)
                            divider = make_axes_locatable(ax)
                            cax = divider.append_axes("right", size = "5%", pad = 0.05)
                            plt.colorbar(im, cax = cax)
                            fname = 'B{}_MOD{}_{}'.format(b, am, 'rho_env_toa_est')
                            plt.savefig('{}/{}_{}.png'.format(output, oname, fname), dpi = 300, bbox_inches = 'tight', facecolor = 'white')
                            plt.close()

                            if (rho_a_opt_mx > rho_a_opt_mn):
                                ## plot bratio
                                fig, ax = plt.subplots()
                                ax.set_title('bratio B{} MOD{}'.format(b, am));
                                im = ax.imshow(bratio, vmin=np.nanpercentile(bratio, 0), vmax=np.nanpercentile(bratio,5))
                                divider = make_axes_locatable(ax)
                                cax = divider.append_axes("right", size = "5%", pad = 0.05)
                                plt.colorbar(im, cax = cax)
                                ax.scatter(minidx[1][0], minidx[0][0], c = 'red', marker = '+')
                                fname = 'B{}_MOD{}_{}'.format(b, am, 'bratio')
                                plt.savefig('{}/{}_{}.png'.format(output, oname, fname), dpi = 300, bbox_inches = 'tight', facecolor = 'white')
                                plt.close()
                        ## END DEVELOPMENT BLOCK ##

                print('End band {}'.format(b))
                print()

            ## START DEVELOPMENT BLOCK ##
            if setu['radcor_development']:
                dev_print['rho_toa_av'] = rho_toa_av
                dev_print['rho_dark'] = rho_dark
                dev_print['psf_tsdsf_mix_C_cov'] = psf_tsdsf_mix_C_cov
                dev_print['psf_tsdsf_mix_M_cov'] = psf_tsdsf_mix_M_cov
                for ai, am in enumerate(aer_models):
                    dev_print['aot550_est_' + am] = tau550_est[ai,]
            ## END DEVELOPMENT BLOCK ##

            #
            # End Run TS-DSF
            # --------------------------------------------------------------

        elif setu['radcor_aot_estimate'] == 'optimise':
            # --------------------------------------------------------------
            # Run optimisation to field measurement:
            # Tests have been done above for number of bands and whether target is in scene
            #

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
                    print('Optimising aot to {} with {} pixel radius and centre pixel: {}, {}'.format(setu['optimise_target_type'], radius_pix, i, j))
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
                    gem.close()
                    #gemo.close()
                    return

            ## cost function for aot estimate
            optimise_aot_cost = setu['optimise_aot_cost'].upper()
            print('Optimising aot with cost function: {}'.format(optimise_aot_cost))

            ## aot optimisation function
            ## QV 2024-05-21 added MAPD
            ## QV 2024-05-22 new version using temporary gem
            def opt_aot(aot, return_res = False, correct_band = correct_band):
                print('    {} {}'.format(lut, aot))
                ## make temporary gem
                ## is not really needed if no GC requested - but
                gemt = ac.gem.gem(None, new = False)
                gemt.gatts = {k: gem.gatts[k] for k in gem.gatts}
                gemt.gatts['uoz'] = uoz
                gemt.gatts['uwv'] = uwv
                gemt.gatts['pressure'] = pressure
                gemt.gatts['wind'] = wind

                ## run through bands doing RAdCor for current model and aot
                opt_bands_ = [b for b in opt_bands]
                if setu['dsf_residual_glint_correction']: opt_bands_ = [b for b in bands_] ##
                for bi, b in enumerate(opt_bands_):
                    ds = bands[b]['rhos_ds']
                    rhos_ = correct_band(b, aot, new = False, write = False)
                    #print(b, ds, rhos_.shape)
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
                    #del gemt.data_mem[ds] ## remove temporary data
                del gemt ## delete temporary gem

                if return_res:
                    return(res_rhos)
                else:
                    if optimise_aot_cost == 'RMSD':
                        return(np.nanmean((opt_rhos-res_rhos)**2)**0.5)
                    elif optimise_aot_cost == 'MAPD':
                        #mn = 0.5 * (opt_rhos + res_rhos) ## MARD
                        mn = 1.0 * opt_rhos ## MAPD
                        mn2 = np.abs(opt_rhos-res_rhos) / mn
                        return(np.nansum(mn2) / len(opt_bands))

            # ## offsets to extract centre of the image
            # cen_offset_0 = int(id_psf)
            # cen_offset_1 = int(id_psf)
            # ## change offsets if radcor_aot_estimate_centre_extent is given
            # if setu['radcor_aot_estimate_centre_extent'] is not None:
            #     cen_offset_npix = np.round(float(setu['radcor_aot_estimate_centre_extent'] * 1000) / resolution).astype(int)
            #     cen_offset_0 = np.round((x_a_dim[0] - cen_offset_npix) / 2).astype(int)
            #     cen_offset_1 = np.round((x_a_dim[1] - cen_offset_npix) / 2).astype(int)
            #     print('Limiting TS-DSF scene centre to radcor_aot_estimate_centre_extent={:.1f} km ({} pixels)'.format(setu['radcor_aot_estimate_centre_extent'], cen_offset_npix))
            # else:
            #     print('Limiting TS-DSF scene to PSF coverage with radius of {:.1f} km ({} pixels)'.format(setu['radcor_kernel_radius'], psf_radius_pixels))
            # print('Offsets to extract image centre: {}x{} pixels'.format(cen_offset_0, cen_offset_1))
            # print('Extracted centre size: {}x{} pixels'.format(x_a_dim[0] - 2*cen_offset_0, x_a_dim[1] - 2 * cen_offset_1))
            ## end offsets to extract centre of the image

            ## Create x vector for RAdCor
            xvec = np.arange(resolution / 2000.0, setu['radcor_kernel_radius'], step = resolution / 1000.0)  # AC 20240306
            xvec = np.hstack((-1.0 * np.flip(xvec), xvec))[1:]                                            # AC 20240306
            id_psf = int((len(xvec) - 1) / 2)                                                             # AC 20240306
            x_f_dim = x_a_dim * 1                                                                         # AC 20240306
            x_o_dim = x_a_dim * 1                                                                         # AC 20240306
            if not setu['radcor_edge_extend']:                                                                 # AC 20240306
                x_o_dim = int(x_a_dim[0] - id_psf * 2), int(x_a_dim[1] - id_psf * 2)                      # AC 20240306
            else:                                                                                         # AC 20240306
                x_f_dim = int(x_a_dim[0] + id_psf * 2), int(x_a_dim[1] + id_psf * 2)                      # AC 20240306

            ## Create uniform PSF and OTF
            if (setu['radcor_kernel_complete_method'] == "neighbourhood"): # AC 20240306
                psf_uni = np.zeros((len(xvec), len(xvec))) + 1.0 / len(xvec)**2                           # AC 20240306
                #otf_uni = np.zeros(x_f_dim)                                                               # AC 20240306 ## AC 2024-09-09 ## psf_uni now added to psf_mix
                #otf_uni[0:psf_uni.shape[0], 0:psf_uni.shape[1]] = psf_uni                                 # AC 20240306 ## AC 2024-09-09 ##
                #otf_uni = np.fft.fftn(otf_uni)                                                            # AC 20240306 ## AC 2024-09-09 ##

            ## Fit AOT(550) for each model
            print('    Estimating aot550 for models {}'.format(', '.join(aer_models)))
            model_band_selection = {}
            for ai, am in enumerate(aer_models):
                print('    Estimating aot550 for model {}'.format(am))
                lut = [lut for lut in luts if aer_nm[am] in lut][0]
                opt = scipy.optimize.minimize_scalar(opt_aot, bounds = (0.001, 5.0), method = 'bounded', options = {"xatol":setu['optimise_tolerance']}) ## AC 2024-09-13 ##
                model_band_selection[am] = {'aot': opt.x, 'fit': opt_aot(opt.x), 'result': opt_aot(opt.x, return_res = True)}
                print('    Optimised aot550 for model {}: {:.4f}'.format(am, model_band_selection[am]['aot']))
                print('    Fit: {:.4f}'.format(model_band_selection[am]['fit']))

            ## plot results
            if setu['optimise_plot']:
                opt_wave = np.asarray([bands[b]['wavelength'] for b in opt_bands])
                opt_sort = np.argsort(opt_wave)

                ## sel model for annotating plot
                sel_am = None
                for am in model_band_selection:
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
                    if am == 'M':
                        col = '#1f77b4'
                    if am == 'C':
                        col = '#ff7f0e'
                    if am == 'U': ## AC 2024-08-27 ##
                        col = '#ff0000' ## AC 2024-08-27 ##
                    plt.plot(opt_wave[opt_sort], model_band_selection[am]['result'][opt_sort],  '.:', color = col,
                            label = r'MOD{} $\tau_a$={:.3f} {}={:.2e} {}'.format(am,
                                                                              model_band_selection[am]['aot'],
                                                                              optimise_aot_cost,
                                                                              model_band_selection[am]['fit'],
                                                                              '(*)' if sel_am == am else ''))
                plt.legend()
                plt.title('{} {} RAdCor'.format(sensor.replace('_', '/'), gem.gatts['isodate'][0:19]))
                plt.xlabel('Wavelength (nm)')
                plt.ylabel(r'$\rho_{s}$ (1)')
                xlim = plt.xlim()
                plt.plot(xlim, [0,0], '--', color='Grey')
                plt.xlim(xlim)
                plt.ylim(setu['optimise_plot_range'][0], setu['optimise_plot_range'][1])
                fname = 'rhos_optimised'
                plt.savefig('{}/{}_{}.png'.format(output, oname, fname), dpi = 300, bbox_inches = 'tight', facecolor = 'white')
                plt.close()

            #
            # End run optimisation to field measurement
            # --------------------------------------------------------------


        #
        # End AOT estimation
        # --------------------------------------------------------------

        # --------------------------------------------------------------
        # Select aerosol model and AOT(550):
        #

        ## fit and aot are already computed for radcor_aot_estimate
        if (setu['radcor_aot_estimate'] == 'optimise') &  (setu['radcor_force_aot'] is None):
            print('\nSelecting best fitting model from optimisation')
            best_mod = None
            best_band = np.nan
            for ai, am in enumerate(aer_models):
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

        else:
            ## track model information
            model_band_selection = {}

            ### select aot per model
            print('Selecting aot per model using radcor_aot_selection={}'.format(setu['radcor_aot_selection']))

            ## select lowest aot
            if setu['radcor_aot_selection'] == 'minimum':
                for ai, am in enumerate(aer_models):
                     lut = [lut for lut in luts if aer_nm[am] in lut][0]
                     min_pos_idx = np.argsort(tau550_est[ai, :][np.where(np.isfinite(tau550_est[ai, :]))])[0]
                     model_band_selection[am] = {'aot': tau550_est[ai, min_pos_idx], 'band': bands_[min_pos_idx],
                                                  'min_pos_idx': min_pos_idx}
                     print('    Model {} estimated AOT(550): {:.4f}'.format(am, model_band_selection[am]['aot']))

            ## select highest aot with lowest positive difference to rho dark
            elif setu['radcor_aot_selection'] == 'lowest_positive_difference':
                for ai, am in enumerate(aer_models):
                    lut = [lut for lut in luts if aer_nm[am] in lut][0]

                    diffs = []

                    for bi,b in enumerate(bands_):
                        if add_rsky:
                            romix_band = np.asarray([lutdw[lut]['rgi'][b_]((pressure, lutdw[lut]['ipd'][romix_par], raa, vza, sza, wind, tau550_est[ai,bi])) for b_ in bands_])
                            # use pure romix
                            #romix_band = np.asarray([lutdw[lut]['rgi'][b_]((pressure, lutdw[lut]['ipd']['romix'], raa, vza, sza, wind, tau550_est[ai,bi])) for b_ in bands_])
                        else:
                            romix_band = np.asarray([lutdw[lut]['rgi'][b_]((pressure, lutdw[lut]['ipd']['romix'], raa, vza, sza, tau550_est[ai,bi])) for b_ in bands_])
                        #min_diff = np.nanmin(rho_dark - romix_band)      #AC_20240221
                        min_diff = np.nanmin(rho_a_est[ai,] - romix_band) #AC_20240221
                        diffs.append(min_diff)

                    diffs = np.asarray(diffs)
                    diffs[diffs<0] = np.nan
                    min_pos_diff = np.nanmin(diffs)
                    if np.isnan(min_pos_diff): ## select minimum if all are negative
                        min_pos_idx = np.argsort(tau550_est[ai,:])[0]
                    else:
                        min_pos_idx = np.where(diffs == min_pos_diff)[0][0]
                    print('Band index, aot', bands_[min_pos_idx], tau550_est[ai,min_pos_idx])

                    model_band_selection[am] = {'aot': tau550_est[ai,min_pos_idx], 'band': bands_[min_pos_idx],
                                                'min_pos_diff': min_pos_diff, 'min_pos_idx': min_pos_idx}

            ### select best fitting model based on two best fitting bands rho_dark:rho_path
            print('Selecting best fitting model')

            ## select model and aot
            if (setu['radcor_force_model'] is not None): ## user supplied model
                best_mod = setu['radcor_force_model']
                best_fit = np.nan
                best_idx = [ai for ai, amod in enumerate(aer_models) if best_mod == amod]
                best_band = np.nan
                print('\nForced aerosol model: {}'.format(best_mod))
                ## check if AOT was also forced
                if setu['radcor_force_aot'] is not None: ## user supplied aot
                    best_aot = 1.0 * setu['radcor_force_aot']
                    print('Forced aerosol optical thickness at 550 nm: {:.4f}'.format(best_aot))
                else: ## aot determined above
                    best_aot = model_band_selection[best_mod]['aot']
                    print('Estimated aerosol optical thickness at 550 nm: {:.4f}'.format(best_aot))
            else: ## model selection
                best_mod = None
                for ai, am in enumerate(aer_models):
                    print('    Model {} estimated AOT(550): {:.4f} (band {})'.format(am, model_band_selection[am]['aot'], model_band_selection[am]['band']))
                    lut = [lut for lut in luts if aer_nm[am] in lut][0]

                    if add_rsky:
                        rho_a_pred = np.asarray([lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd'][romix_par],\
                                                                    raa, vza, sza, wind, model_band_selection[am]['aot'])) for i, b in enumerate(bands_)])
                        # use pure romix
                        #rho_a_pred = np.asarray([lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd'][romix_par],\
                        #                                            raa, vza, sza, wind, model_band_selection[am]['aot'])) for i, b in enumerate(bands_)])
                    else:
                        rho_a_pred = np.asarray([lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['romix'],\
                                                                    raa, vza, sza, model_band_selection[am]['aot'])) for i, b in enumerate(bands_)])

                    ## full spectrum
                    #model_fit = (np.nanmean(abs(rho_dark - rho_a_pred) / rho_dark))

                    ## best two bands
                    #diff = abs(rho_dark - rho_a_pred) / rho_dark       # AC 20240221
                    diff = abs(rho_a_est[ai,] - rho_a_pred) / rho_dark  # AC 20240221
                    bidx = np.argsort(diff)
                    #print(diff)
                    # print(bands_[bidx[0]], diff[bidx[0]])
                    # print(bands_[bidx[1]], diff[bidx[1]])
                    # model_fit = (np.nanmean((abs(rho_dark - rho_a_pred) / rho_dark)[bidx[0:2]]))             # AC 20240221
#                    model_fit = (np.nanmean((abs(rho_a_est[ai,] - rho_a_pred) / rho_a_est[ai,])[bidx[0:2]]))   # AC 20240221
                    tau_eval = 1 * tau550_est # ALEX 2024-08-27
                    tau_eval[tau_eval == 5] = np.nan # ALEX 2024-08-27
                    model_fit = np.nanstd( tau_eval[ai,] ) / np.nanmean( tau_eval[ai,] ) # ALEX 2024-08-27
#                    model_fit = np.nansum( tau_eval[ai,] / np.nanmin(tau_eval[ai,]) ) # ALEX 2024-08-27

#                    print('    Model {} rho_a fit to rho_a_est: {:.4f}'.format(am, model_fit)) # ALEX 2024-08-27
                    print('    Model {} tau_aer_550 CV: {:.4f}'.format(am, model_fit)) # ALEX 2024-08-27

                    if best_mod is None:
                        best_mod = am
                        best_fit = 1.0 * model_fit
                        best_idx = ai
                        best_band = bidx[0]
                        best_aot = model_band_selection[am]['aot']
                    elif best_fit > model_fit:
                        best_mod = am
                        best_fit = 1.0 * model_fit
                        best_idx = ai
                        best_band = bidx[0]
                        best_aot = model_band_selection[am]['aot']

                print('\nEstimated aerosol model: {} (fitted band {})'.format(best_mod, bands_[best_band]))
                print('Estimated aerosol optical thickness at 550 nm: {:.4f}'.format(best_aot))

            ## START DEVELOPMENT BLOCK ##
            if setu['radcor_development']:
                dev_print['selected_model'] = best_mod
                dev_print['selected_aot550'] = best_aot
            ## END START DEVELOPMENT BLOCK ##

            if setu['radcor_diagnostic_plots']:
                ## plot difference between estimated and modeled path reflectance
                for ai, am in enumerate(aer_models):
                    fig, ax = plt.subplots()
                    lut = [lut for lut in luts if aer_nm[am] in lut][0]
                    best_model_band = None
                    best_model_band_diff = 0
                    for bi,b in enumerate(bands_):
                        if add_rsky:
                            romix_band = np.asarray([lutdw[lut]['rgi'][b_]((pressure, lutdw[lut]['ipd'][romix_par], raa, vza, sza, wind, tau550_est[ai,bi])) for b_ in bands_])
                            ## use pure_romix
                            #romix_band = np.asarray([lutdw[lut]['rgi'][b_]((pressure, lutdw[lut]['ipd']['romix'], raa, vza, sza, wind, tau550_est[ai,bi])) for b_ in bands_])
                        else:
                            romix_band = np.asarray([lutdw[lut]['rgi'][b_]((pressure, lutdw[lut]['ipd']['romix'], raa, vza, sza, tau550_est[ai,bi])) for b_ in bands_])
                        min_diff = np.nanmin(rho_a_est[ai,] - romix_band)
                        alpha = 1 if (best_idx == ai) & (best_band == bi) else 0.25
                        plt.plot(rho_a_est[ai,] - romix_band, '--', label = '{} {:.2f}'.format(b, tau550_est[ai,bi]), alpha=alpha)
                    plt.legend(title=r'$\tau_a$ $550~nm$')
                    plt.title(am)
                    xlim = plt.xlim()
                    plt.plot(xlim, [0,0], '-', color='Black', zorder = 0)
                    plt.xlim(xlim)
                    plt.xlabel('Band index (1)')
                    plt.ylabel(r'$\rho_{path}$ estimated-modeled (1)')
                    fname = 'MOD{}_rho_path_difference'.format(am)
                    plt.savefig('{}/{}_{}.png'.format(output, oname, fname), dpi = 300, bbox_inches = 'tight', facecolor = 'white')
                    plt.close()


                ## plot selected modeled path reflectance and estimated path reflectance
                for ai, am in enumerate(aer_models):
                    fig, ax = plt.subplots()
                    lut = [lut for lut in luts if aer_nm[am] in lut][0]
                    for bi,b in enumerate(bands_):
                        if add_rsky:
                            romix_band = [lutdw[lut]['rgi'][b_]((pressure, lutdw[lut]['ipd'][romix_par], raa, vza, sza, wind, tau550_est[ai,bi])) for b_ in bands_]
                            # use pure romix
                            #romix_band = [lutdw[lut]['rgi'][b_]((pressure, lutdw[lut]['ipd']['romix'], raa, vza, sza, wind, tau550_est[ai,bi])) for b_ in bands_]
                        else:
                            romix_band = [lutdw[lut]['rgi'][b_]((pressure, lutdw[lut]['ipd']['romix'], raa, vza, sza, tau550_est[ai,bi])) for b_ in bands_]
                        alpha = 1 if (best_idx == ai) & (best_band == bi) else 0.25
                        plt.plot(romix_band, '--', alpha=alpha, label = '{} {:.2f}'.format(b, tau550_est[ai,bi]))
                    plt.legend(title=r'$\tau_a$ $550~nm$')
                    plt.title(am)
                    plt.plot(rho_a_est[ai,], color='Black')
                    plt.xlabel('Band index (1)')
                    plt.ylabel(r'$\rho_{path}$ (1)')
                    fname = 'MOD{}_rho_path_selected'.format(am)
                    plt.savefig('{}/{}_{}.png'.format(output, oname, fname), dpi = 300, bbox_inches = 'tight', facecolor = 'white')
                    plt.close()

                ## plot taua spectra estimated per band
                fig, ax = plt.subplots()
                for ai, am in enumerate(aer_models):
                    if am == 'M':
                        col = '#1f77b4'
                    if am == 'C':
                        col = '#ff7f0e'
                    if am == 'U':
                        col = '#ff0000'
                    plt.plot(tau550_est[ai,:], label = 'MOD{}'.format(am), color = col)
                plt.legend()
                plt.xlabel('Band index (1)')
                plt.ylabel(r'$\tau_{a}$ $550~nm$ (1)')
                fname = 'taua_bands'
                plt.savefig('{}/{}_{}.png'.format(output, oname, fname), dpi = 300, bbox_inches = 'tight', facecolor = 'white')
                plt.close()

    #
    # End Select aerosol model and AOT(550)
    # --------------------------------------------------------------

    # --------------------------------------------------------------
    # Atmospheric correction proper:
    #

    nbands_ = len(bands_)

    ## START DEVELOPMENT BLOCK ##
    if setu['radcor_development']:
        saf_mix_cov = np.zeros(nbands_)
        psf_mix_cov = np.zeros(nbands_)
    ## END DEVELOPMENT BLOCK ##

    print('\nPerforming RAdCor correction')
    #am = '{}'.format(aer_models[best_idx]) # AC 20240305
    am = best_mod                           # AC 20240305
    tau_aer_550 = best_aot * 1.0
    lut = [lut for lut in luts if aer_nm[am] in lut][0]

    print('taer_550: {}'.format(tau_aer_550))
    print('Aerosol model: {}'.format(am))
    print('Aerosol LUT: {}'.format(lut))

    ## Create x vector for RAdCor
    xvec = np.arange(resolution / 2000.0, setu['radcor_kernel_radius'], step = resolution / 1000.0)  # AC 20240306
    xvec = np.hstack((-1.0 * np.flip(xvec), xvec))[1:]                                            # AC 20240306
    id_psf = int((len(xvec) - 1) / 2)                                                             # AC 20240306
    x_f_dim = x_a_dim * 1                                                                         # AC 20240306
    x_o_dim = x_a_dim * 1                                                                         # AC 20240306
    if not setu['radcor_edge_extend']:                                                                 # AC 20240306
        x_o_dim = int(x_a_dim[0] - id_psf * 2), int(x_a_dim[1] - id_psf * 2)                      # AC 20240306
    else:                                                                                         # AC 20240306
        x_f_dim = int(x_a_dim[0] + id_psf * 2), int(x_a_dim[1] + id_psf * 2)                      # AC 20240306

    # ## Create scene edge mask based on RAdCor PSF
    scene_edge_mask = np.zeros((x_a_dim[0], x_a_dim[1])) + np.nan
    scene_edge_mask[id_psf:x_a_dim[0]-id_psf, id_psf:x_a_dim[1]-id_psf] = 1.0

    ## 2024-12-02
    ## offsets to extract centre of the image
    cen_offset_0 = int(id_psf)
    cen_offset_1 = int(id_psf)
    ## change offsets if radcor_aot_estimate_centre_extent is given
    if setu['radcor_aot_estimate_centre_extent'] is not None:
        cen_offset_npix = np.round(float(setu['radcor_aot_estimate_centre_extent'] * 1000) / resolution).astype(int)
        cen_offset_0 = np.round((x_a_dim[0] - cen_offset_npix) / 2).astype(int)
        cen_offset_1 = np.round((x_a_dim[1] - cen_offset_npix) / 2).astype(int)
    ## 2024-12-02

    ## Create uniform PSF and OTF
    if (setu['radcor_kernel_complete_method'] == "neighbourhood"): # AC 20240306
        psf_uni = np.zeros((len(xvec), len(xvec))) + 1.0 / len(xvec)**2                           # AC 20240306
        #otf_uni = np.zeros(x_f_dim)                                                               # AC 20240306 ## AC 2024-09-09 ## psf_uni now added to psf_mix
        #otf_uni[0:psf_uni.shape[0], 0:psf_uni.shape[1]] = psf_uni                                 # AC 20240306
        #otf_uni = np.fft.fftn(otf_uni)                                                            # AC 20240306

    # --------------------------------------------------------------
    # Create ACOLITE output file:
    #

    #if setu['radcor_data_in_memory']: gem.store = True ## not faster in limited tests

    ## set up output gem
    gemo = ac.gem.gem(ofile, new = True)
    gemo.gatts = {k: gem.gatts[k] for k in gem.gatts}
    gemo.gatts['uoz'] = uoz
    gemo.gatts['uwv'] = uwv
    gemo.gatts['pressure'] = pressure
    gemo.gatts['wind'] = wind
    gemo.gatts['ofile'] = ofile

    ## Add rhoe to BEAM format auto-grouping
    gemo.gatts['auto_grouping'] = 'rhot:rhorc:rhos:rhow:Rrs:rhoe:rhosu'
    if setu['radcor_write_rhotc']: gemo.gatts['auto_grouping']+=':rhotc'
    if setu['output_ed']: gemo.gatts['auto_grouping']+=':Ed' # :Edu:Ed0
    gemo.gatts['acolite_file_type'] = 'L2R'
    gemo.gatts['radcor_version'] = '{}'.format(ac.adjacency.radcor.version)
    if gem.nc_projection is not None:
        gemo.nc_projection = gem.nc_projection
        if setu['radcor_crop_centre']: ## crop to centre area
            gemo.nc_projection['x']['data'] = gemo.nc_projection['x']['data'][cen_offset_1:x_a_dim[1] - cen_offset_1]
            gemo.nc_projection['y']['data'] = gemo.nc_projection['y']['data'][cen_offset_0:x_a_dim[0] - cen_offset_0]
    gemo.bands = bands

    ## add selected model and aot
    gemo.gatts['ac_model'] = best_mod
    gemo.gatts['ac_aot_550'] = best_aot
    gemo.gatts['ac_lut'] = lut
    gemo.gatts['ac_fit'] = best_fit

    ## store bands used for TS-DSF
    gemo.gatts['ac_bands'] = ','.join([str(b) for b in bands_])
    ## store selected band index and name (only one band since we choose the model spectrally)
    if np.isfinite(best_band):
        gemo.gatts['ac_band1_idx'] = best_band
        gemo.gatts['ac_band1'] = bands_[best_band]

    ## also create separate L1R file for rhotc if radcor_write_rhotc_separate_file
    if setu['radcor_write_rhotc'] & setu['radcor_write_rhotc_separate_file']:
        ofile_l1rc = '{}/{}.nc'.format(output, oname.replace('L2R', 'L1RC'))
        gemo_l1rc = ac.gem.gem(ofile_l1rc, new = True)
        ## copy gatts
        gemo_l1rc.gatts = {k: gemo.gatts[k] for k in gemo.gatts}
        gemo_l1rc.gatts['acolite_file_type'] = 'L1RC'
        ## copy projection and bands
        gemo_l1rc.nc_projection = gemo.nc_projection
        gemo_l1rc.bands = gemo.bands

    ## Read and store datasets - add cropping step
    sub = None
    for ds in gem.datasets:
        if ds not in setu['copy_datasets']: continue
        if 'projection_key' in gem.gatts:
            if ds in ['x', 'y', gem.gatts['projection_key']]: continue
        print('Reading {} from inputfile {}'.format(ds, ncf))
        d_, a_ = gem.data(ds, attributes = True)
        print('Writing {} to outputfile {}'.format(ds, ofile))
        if setu['radcor_crop_centre']: ## crop to centre area
            d_ = d_[cen_offset_0:x_a_dim[0] - cen_offset_0, cen_offset_1:x_a_dim[1] - cen_offset_1]
        gemo.write(ds, d_, ds_att = a_)
        ## write also to L1RC file
        if setu['radcor_write_rhotc'] & setu['radcor_write_rhotc_separate_file']:
            gemo_l1rc.write(ds, d_, ds_att = a_)

    if setu['dem_pressure'] & setu['dem_pressure_write'] & (elevation is not None):
        if setu['radcor_crop_centre']: ## crop to centre area
            elevation = elevation[cen_offset_0:x_a_dim[0] - cen_offset_0, cen_offset_1:x_a_dim[1] - cen_offset_1]
        gemo.write('dem', elevation)
        gemo.write('dem_pressure',ac.ac.pressure_elevation(elevation))
        elevation = None
        gemo.gatts['elevation'] = median_elevation

    #
    # End Create ACOLITE output file
    # --------------------------------------------------------------


    ## Run through bands
    for bi, b in enumerate(bands):
        if not bands[b]['radcor_use_band']:
            if setu['radcor_write_rhot']:
                print('Loading TOA data for band {} ({})'.format(b, bands[b]['rhot_ds']))
                rho_toa, att = gem.data(bands[b]['rhot_ds'], attributes = True)
                print('Writing band {}'.format(b))
                if setu['radcor_crop_centre']:
                    rho_toa = rho_toa[cen_offset_0:x_a_dim[0] - cen_offset_0, cen_offset_1:x_a_dim[1] - cen_offset_1]
                gemo.write(bands[b]['rhot_ds'], rho_toa, ds_att = att)
                ## write to L1RC file as well
                if setu['radcor_write_rhotc'] & setu['radcor_write_rhotc_separate_file']:
                    gemo_l1rc.write(bands[b]['rhot_ds'], rho_toa, ds_att = att)
                del rho_toa
                del att
            continue
        ## run correct band function
        correct_band(b, tau_aer_550, write = True, quiet = False)
    ## End bands loop
    print('Wrote {}'.format(ofile))

    ## START DEVELOPMENT BLOCK ##
    if setu['radcor_development']:
        print('(DEV) PSF radius: {} (pixels)'.format(psf_radius_pixels))
        print('(DEV) PSF dimension: {}x{} (pixels)'.format(int(psf_radius_pixels * 2), int(psf_radius_pixels * 2)))
        print('(DEV) Extended PSF dimensions: {}x{} (pixels)'.format(x_f_dim[0], x_f_dim[1]))
        print('(DEV) Length of xvec: {}'.format(len(xvec)))
        print('(DEV) Index id_psf: {}'.format(id_psf))
        print('\n(DEV) Breakpoints for the PSF calculation: \n{}'.format(xvec))
        dev_print['psf_pixel_radius'] = psf_radius_pixels
        dev_print['psf_pixel_dim'] = int(psf_radius_pixels * 2), int(psf_radius_pixels * 2)
        dev_print['id_psf'] = id_psf
        dev_print['xvec_len'] = len(xvec)
        dev_print['xvec'] = xvec

        dev_print['saf_mix_cov'] = saf_mix_cov
        dev_print['psf_mix_cov'] = psf_mix_cov

        ## Get atmospheric parameters
        rho_a     = np.zeros(nbands_)
        rho_a_sph = np.zeros(nbands_)
        T_d_tot   = np.zeros(nbands_)
        T_u_tot   = np.zeros(nbands_)
        T_u_dir   = np.zeros(nbands_)
        T_u_dif   = np.zeros(nbands_)
        tau_tot   = np.zeros(nbands_)

        ## Run through bands
        for bi, b in enumerate(bands_):
            if add_rsky:
                rho_a[bi]     = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd'][romix_par], raa, vza, sza, wind, tau_aer_550))
                rho_a_sph[bi] = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['astot'], raa, vza, sza, wind, tau_aer_550))
                T_d_tot[bi]   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['dtott'], raa, vza, sza, wind, tau_aer_550))
                T_u_tot[bi]   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['utott'], raa, vza, sza, wind, tau_aer_550))
                tau_tot[bi]   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['ttot'],  raa, vza, sza, wind, tau_aer_550))
            else:
                rho_a[bi]     = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['romix'], raa, vza, sza, tau_aer_550))
                rho_a_sph[bi] = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['astot'], raa, vza, sza, tau_aer_550))
                T_d_tot[bi]   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['dtott'], raa, vza, sza, tau_aer_550))
                T_u_tot[bi]   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['utott'], raa, vza, sza, tau_aer_550))
                tau_tot[bi]   = lutdw[lut]['rgi'][b]((pressure, lutdw[lut]['ipd']['ttot'],  raa, vza, sza, tau_aer_550))
            T_u_dir[bi]   = np.exp(-tau_tot[bi] / cos_vza)
            T_u_dif[bi]   = T_u_tot[bi] - T_u_dir[bi]

        dev_print['tau_tot'] = tau_tot
        dev_print['rho_a'] = rho_a
        dev_print['rho_a_sph'] = rho_a_sph
        dev_print['T_d_tot'] = T_d_tot
        dev_print['T_u_tot'] = T_u_tot
        dev_print['T_u_dif'] = T_u_dif
        dev_print['T_u_dir'] = T_u_dir

    ## END DEVELOPMENT BLOCK ##

    #
    # End Atmospheric correction proper
    # --------------------------------------------------------------

    ## update dataset info
    gemo.setup()

    ## glint correction
    if setu['dsf_residual_glint_correction']:
        if setu['dsf_residual_glint_correction_method'] != 'default':
            print('dsf_residual_glint_correction_method={} not implemented after RAdCor'.format(setu['dsf_residual_glint_correction_method']))
        else:
            print('Running glint correction!')
            ## run glint correction
            ret = ac.glint.default(gemo, settings = setu, new_file = False, write = True, lutdw = lutdw)

    ## close files
    gem, gemo = None, None
    if setu['radcor_development']:
        gempsf, gemconf = None, None
    gemo_l1rc = None

    return(ofile)
