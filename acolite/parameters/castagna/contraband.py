## def contraband
## Computes Castagna et al. 2020 orange contraband  https://doi.org/10.3390/rs12040637
## and other contrabands
##
## written by Quinten Vanhellemont, RBINS
## 2022-06-21
## modifications: 2023-03-05 (QV) renamed from orange and added other contrabands for L7 and Pléiades
##

def contraband(gem, verbosity=5):
    import acolite as ac
    import os

    ## read gem file if NetCDF
    if type(gem) is str:
        gem = ac.gem.gem(gem)
        nc_projection = gem.nc_projection
    gemf = gem.file

    ## compute contrabands
    cb_file = ac.config['data_dir']+'/Shared/algorithms/Castagna/{}_contra_coefficients.txt'.format(gem.gatts['sensor'])
    if os.path.exists(cb_file):
        ## load config for this sensor
        cb_cfg = ac.shared.import_config(cb_file)
        cb_cfg = {k:cb_cfg[k].split(',') for k in cb_cfg}
        print('Computing contrabands: {}'.format(', '.join(cb_cfg['contrabands'])))

        ## required bands
        req_bands = cb_cfg['pan'] + cb_cfg['ms']
        sensor_cb = '{}_CONTRA'.format(gem.gatts['sensor'])

        ## make bands dataset
        rsrd = ac.shared.rsr_dict(gem.gatts['sensor'])[gem.gatts['sensor']]
        gem.bands = {}
        for bi, b in enumerate(rsrd['rsr_bands']):
            if b not in gem.bands:
                gem.bands[b] = {k:rsrd[k][b] for k in ['wave_mu', 'wave_nm', 'wave_name'] if b in rsrd[k]}
                gem.bands[b]['rhot_ds'] = 'rhot_{}'.format(gem.bands[b]['wave_name'])
                gem.bands[b]['rhos_ds'] = 'rhos_{}'.format(gem.bands[b]['wave_name'])
                gem.bands[b]['wavelength']=gem.bands[b]['wave_nm']
        ## end bands dataset

        ## do we have the required bands
        datasets = []
        compute_contraband=True
        for b in req_bands:
            if gem.bands[b]['rhos_ds'] not in gem.datasets:
                print('{} not present, skipping contraband computation'.format(gem.bands[b]['rhos_ds']))
                compute_contraband=False
            else:
                datasets.append(gem.bands[b]['rhos_ds'])

        if compute_contraband:
            ## read rsr for wavelength name
            rsrd_cb = ac.shared.rsr_dict(sensor_cb)[sensor_cb]

            for b in cb_cfg['contrabands']:
                ## set up band attributes
                cb = {k:rsrd_cb[k][b] for k in ['wave_mu', 'wave_nm', 'wave_name']}
                cb['rhos_ds'] = 'rhos_{}'.format(cb['wave_name'])
                cb['wavelength'] = cb['wave_nm']
                gem.bands[b] = cb
                print('Computing contraband {}'.format(cb['rhos_ds']))

                ## compute contraband
                for ii in range(len(datasets)):
                    if ii == 0:
                        cb_data = gem.data(datasets[ii])*float(cb_cfg['pf_{}'.format(b.lower())][0])
                    else:
                        cb_data += gem.data(datasets[ii])*float(cb_cfg['msf_{}'.format(b.lower())][ii-1])
                gem.write(cb['rhos_ds'], cb_data, ds_att = cb)
                cb_data = None
                cb = None
        ## end contraband
    else:
        print('No contrabands configured for {}.'.format(gem.gatts['sensor']))
