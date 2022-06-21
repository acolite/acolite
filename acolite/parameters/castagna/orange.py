## def orange
## Computes Castagna et al. 2020 orange band  https://doi.org/10.3390/rs12040637
## split off from acolite_l2r
##
## written by Quinten Vanhellemont, RBINS
## 2022-06-21
## modifications:
##

def orange(gem, verbosity=5):
    import acolite as ac

    ## read gem file if NetCDF
    if type(gem) is str:
        gem = ac.gem.gem(gem)
        nc_projection = gem.nc_projection
    gemf = gem.file

    ## compute oli orange band
    if (gem.gatts['sensor'] in ['L8_OLI', 'L9_OLI', 'EO1_ALI']):
        if verbosity > 1: print('Computing orange band')

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

        ## load orange band configuration
        if gem.gatts['sensor'] == 'L8_OLI':
            ob_cfg = ac.shared.import_config(ac.config['data_dir']+'/L8/oli_orange.cfg')
            sensor_o = 'L8_OLI_ORANGE'
            panb, greenb, redb = '8', '3', '4'
        if gem.gatts['sensor'] == 'L9_OLI':
            ob_cfg = ac.shared.import_config(ac.config['data_dir']+'/L9/oli_orange.cfg')
            sensor_o = 'L9_OLI_ORANGE'
            panb, greenb, redb = '8', '3', '4'
        if gem.gatts['sensor'] == 'EO1_ALI':
            ob_cfg = ac.shared.import_config(ac.config['data_dir']+'/EO1/ali_orange.cfg')
            sensor_o = 'EO1_ALI_ORANGE'
            panb, greenb, redb = '1', '4', '5'

        ## do we have the required bands
        compute_orange = True
        for b in [panb, greenb, redb]:
            if gem.bands[b]['rhos_ds'] not in gem.datasets:
                print('{} not present, skipping orange band computation'.format(gem.bands[b]['rhos_ds']))
                compute_orange=False

        if compute_orange:
            ## read rsr for wavelength name
            rsrd_o = ac.shared.rsr_dict(sensor_o)[sensor_o]
            ob = {k:rsrd_o[k]['O'] for k in ['wave_mu', 'wave_nm', 'wave_name']}
            ob['rhos_ds'] = 'rhos_{}'.format(ob['wave_name'])
            ob['wavelength'] = ob['wave_nm']
            gem.bands['O'] = ob
            ## compute orange band
            ob_data = gem.data(gem.bands[panb]['rhos_ds'])*float(ob_cfg['pf'])
            ob_data += gem.data(gem.bands[greenb]['rhos_ds'])*float(ob_cfg['gf'])
            ob_data += gem.data(gem.bands[redb]['rhos_ds'])*float(ob_cfg['rf'])
            gem.write(ob['rhos_ds'], ob_data, ds_att = ob)
            #print('Wrote orange band {} to {}'.format(ob['rhos_ds'], gemf))
            ob_data = None
            ob = None
    ## end orange band
