## get sensor, rsr and bands dict for NetCDF file
##
## written by Quinten Vanhellemont, RBINS
## 2024-05-15
## modifications: 2024-05-15 (QV) renamed from gem_bands, added as ac.gem.bands
##                2025-01-20 (QV) added rhosu bands, added rsr versioning

def bands(gem):
    import acolite as ac

    ## check if file path is given
    opened = False
    if type(gem) is str:
        gem = ac.gem.gem(gem)
        opened = True

    ## get rsr
    sensor = gem.gatts['sensor']
    if ac.settings['run']['rsr_version'] is not None: sensor = '{}_{}'.format(sensor, ac.settings['run']['rsr_version'])
    rsrd = ac.shared.rsr_dict(sensor=sensor)[sensor]

    ## get ozone/water vapour from defaults or gatts
    uoz = ac.settings['run']['uoz_default']
    uwv = ac.settings['run']['uwv_default']
    if 'uoz' in gem.gatts: uoz = gem.gatts['uoz']
    if 'uwv' in gem.gatts: uwv = gem.gatts['uwv']

    ## get gas transmittance
    tg_dict = ac.ac.gas_transmittance(gem.gatts['sza'], gem.gatts['vza'], uoz = uoz, uwv = uwv, rsr = rsrd['rsr'])

    ## make bands dataset
    bands = {}
    for bi, b in enumerate(rsrd['rsr_bands']):
        if b not in bands:
            bands[b] = {k:rsrd[k][b] for k in ['wave_mu', 'wave_nm', 'wave_name'] if b in rsrd[k]}
            bands[b]['rhot_ds'] = 'rhot_{}'.format(bands[b]['wave_name'])
            bands[b]['rhos_ds'] = 'rhos_{}'.format(bands[b]['wave_name'])
            bands[b]['rhosu_ds'] = 'rhosu_{}'.format(bands[b]['wave_name'])
            if ac.settings['run']['add_band_name']:
                bands[b]['rhot_ds'] = 'rhot_{}_{}'.format(b, bands[b]['wave_name'])
                bands[b]['rhos_ds'] = 'rhos_{}_{}'.format(b, bands[b]['wave_name'])
                bands[b]['rhosu_ds'] = 'rhosu_{}_{}'.format(b, bands[b]['wave_name'])
            if ac.settings['run']['add_detector_name']:
                dsname = rhot_ds[bi][5:]
                bands[b]['rhot_ds'] = 'rhot_{}'.format(dsname)
                bands[b]['rhos_ds'] = 'rhos_{}'.format(dsname)
                bands[b]['rhosu_ds'] = 'rhosu_{}'.format(dsname)
            for k in tg_dict:
                if k not in ['wave']:
                    bands[b][k] = tg_dict[k][b]
                    if ac.settings['run']['gas_transmittance'] is False: bands[b][k] = 1.0
            bands[b]['wavelength']=bands[b]['wave_nm']
    ## end bands dataset

    if opened:
        gem.close()
        del gem

    return(sensor, rsrd, bands)
