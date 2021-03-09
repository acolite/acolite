## def acolite_l2w
## new L2W parameter computation for L2R generic extracted miniscene
## written by Quinten Vanhellemont, RBINS
## 2021-03-09
## modifications:


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

    ## read gem file if NetCDF
    if type(gem) is str:
        gemf = '{}'.format(gem)
        gem = ac.gem.read(gem, sub=sub, load_data=load_data)
    gemf = gem['gatts']['gemfile']

    ## set up output file
    if target_file is None:
        output_name = os.path.basename(gemf).replace('.nc', '')
        output_name = output_name.replace('_L2R', '_L2W')
        odir = output if output is not None else os.path.dirname(gemf)
        ofile = '{}/{}.nc'.format(odir, output_name)
    else:
        ofile = '{}'.format(target_file)

    ## combine default and user defined settings
    setu = ac.acolite.settings.parse(gem['gatts']['sensor'], settings=settings)

    ## get rhot and rhos wavelengths
    rhot_ds = [ds for ds in gem['datasets'] if 'rhot_' in ds]
    rhot_waves = [int(ds.split('_')[1]) for ds in rhot_ds]
    rhos_ds = [ds for ds in gem['datasets'] if 'rhos_' in ds]
    rhos_waves = [int(ds.split('_')[1]) for ds in rhos_ds]
    if len(rhos_waves) == 0: print('{} is probably not an ACOLITE L2R file.'.format(gemf))

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
    if setu['l2w_mask_smooth']: cur_data = scipy.ndimage.gaussian_filter(cur_data, setu['l2w_mask_smooth_sigma'])
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
        if setu['l2w_mask_smooth']: cur_data = scipy.ndimage.gaussian_filter(cur_data, setu['l2w_mask_smooth_sigma'])
        cirrus_mask = cur_data > setu['l2w_mask_cirrus_threshold']
        cirrus = None
        l2_flags += cirrus_mask.astype(np.int32)*(2**setu['flag_exponent_cirrus'])
        cirrus_mask = None
    else:
        if verbosity > 1: print('No suitable band found for cirrus masking.')
    ## TOA out of limit
    toa_mask = None
    for ci, cur_par in enumerate(rhot_ds):
        if cur_par in gem['data']:
            cur_data = 1.0 * gem['data'][cur_par]
        else:
            cur_data = ac.shared.nc_data(gemf, cur_par, sub=sub).data
        if setu['l2w_mask_smooth']: cur_data = scipy.ndimage.gaussian_filter(cur_data, setu['l2w_mask_smooth_sigma'])
        if toa_mask is None: toa_mask = np.zeros(cur_data.shape).astype(bool)
        toa_mask = (toa_mask) | (cur_data > setu['l2w_mask_high_toa_threshold'])
    l2_flags = (l2_flags) | (toa_mask.astype(np.int32)*(2**setu['flag_exponent_toa']))
    toa_mask = None
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
        if (cur_par in setu['l2w_parameters']):
            copy_datasets.append(cur_par)
        if (('rhot_*' in setu['l2w_parameters']) & ('rhot_' in cur_par)):
            copy_datasets.append(cur_par)
        if (('rhos_*' in setu['l2w_parameters']) & ('rhos_' in cur_par)):
            copy_datasets.append(cur_par)
        if (('rhow_*' in setu['l2w_parameters']) & ('rhos_' in cur_par)):
            copy_datasets.append(cur_par.replace('rhos_', 'rhow_'))
        if (('Rrs_*' in setu['l2w_parameters']) & ('rhos_' in cur_par)):
            copy_datasets.append(cur_par.replace('rhos_', 'Rrs_'))

    ## copy datasets
    for ci, cur_par in enumerate(copy_datasets):
        factor = 1.0
        cur_tag = '{}'.format(cur_par)
        mask = False
        ## copy Rrs/rhow
        if 'Rrs_' in cur_par:
            factor = 1.0/np.pi
            cur_tag = cur_par.replace('Rrs_','rhos_')
            mask = True
        if 'rhow_' in cur_par:
            factor = 1.0
            cur_tag = cur_par.replace('rhow_','rhos_')
            mask = True
        ## if data already read copy here
        if cur_tag in gem['data']:
            cur_data = factor * gem['data'][cur_tag]
            cur_att = gem['atts'][cur_tag]
        else:
            cur_d, cur_att = ac.shared.nc_data(gemf, cur_tag, sub=sub, attributes=True)
            cur_data = factor * cur_d.data
            cur_data[cur_d.mask] = np.nan
            cur_d = None
        ## apply mask to Rrs and rhow
        if mask: cur_data[(l2_flags & flag_value)!=0] = np.nan
        if verbosity > 1: print('Writing {}'.format(cur_par))
        ac.output.nc_write(ofile, cur_par, cur_data, dataset_attributes=cur_att, attributes=gem['gatts'], new=new)
        cur_data = None
        new = False

    ## write l2 flags
    ac.output.nc_write(ofile, 'l2_flags', l2_flags, attributes=gem['gatts'], new=new)
    new = False

    ## parameter loop
    ## compute other parameters
    for cur_par in setu['l2w_parameters']:
        if cur_par.lower() in ['rhot_*', 'rhos_*', 'rrs_*', 'rhow_*']: continue ## we have copied these above
        if cur_par.lower() in [ds.lower() for ds in ac.shared.nc_datasets(ofile)]: continue ## parameter already in output dataset (would not work if we are appending subsets to the ncdf)

        ## default mask and empty dicts for current parameter
        mask = False
        par_data = {}
        par_atts = {}

        ## nechad turbidity/spm
        if 'nechad' in cur_par:
            ## turbidity
            if cur_par[0] == 't':
                npar = 'turbidity'
                nechad_dict = ac.parameters.nechad.coef_hyper('TUR')
            elif cur_par[0] == 's':
                npar = 'spm'
                nechad_dict = ac.parameters.nechad.coef_hyper('SPM')
            else:
                continue
        ## end nechad turbidity/spm

        #############################
        ## Pitarch 3 band QAA
        if cur_par[0:5] == 'p3qaa':
            mask = True ## water parameter so apply mask
            ## load config
            cfg = ac.parameters.pitarch.p3qaa_coef()
            if gem['gatts']['sensor'] not in cfg:
                print('P3QAA not configured for {}'.format(gem['gatts']['sensor']))
                continue

            ## read Blue Green Red data, convert to Rrs
            for k in ['B', 'G', 'R']:
                ci, cw = ac.shared.closest_idx(rhos_waves, cfg[gem['gatts']['sensor']]['center_wl'][k])
                cur_ds = 'rhos_{}'.format(cw)
                if cur_ds in gem['data']:
                    cur_data = 1.0 * gem['data'][cur_ds]
                else:
                    cur_data = ac.shared.nc_data(gemf, cur_ds, sub=sub).data
                if mask: cur_data[(l2_flags & flag_value)!=0] = np.nan
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

            ## reformat for output
            for p in cur_par_out:
                if p[6:] not in ret: ## this means we have a three band parameter
                    k = p[6:-4]
                    w = int(p[-3:])
                    if w == Bwave: wi = 0
                    if w == Gwave: wi = 1
                    if w == Rwave: wi = 2
                    par_data[p] = ret[k][:,:,wi]
                    par_atts[p] = {'ds_name':p}
                else:
                    par_data[p] = ret[p[6:]]
                    par_atts[p] = {'ds_name':p}
        ## end Pitarch 3 band QAA
        #############################


        ## continue if parameter not computed
        if len(par_data) == 0:
            print('Parameter {} not computed'.format(cur_par))
            continue

        ## write parameters in par_data
        for cur_ds in par_data:
            ## add mask
            if mask: par_data[cur_ds][(l2_flags & flag_value)!=0] = np.nan
            ## write to NetCDF
            if verbosity > 1: print('Writing {}'.format(cur_ds))
            ac.output.nc_write(ofile, cur_ds, par_data[cur_ds], dataset_attributes=par_atts[cur_ds],
                               attributes=gem['gatts'], new=new)
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
