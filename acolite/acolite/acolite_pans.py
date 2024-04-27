## def acolite_pans
## function to pansharpen ACOLITE NetCDF outputs
## written by Quinten Vanhellemont, RBINS
## 2022-02-20
## modifications: 2022-02-21 (QV) acolite integration
##                2024-03-14 (QV) update settings handling
##                                removed RGB/geotiff outputs
##                2024-04-16 (QV) use new gem NetCDF handling

def acolite_pans(ncf, output = None, settings = None):

    import os
    import acolite as ac
    import numpy as np
    import scipy.ndimage

    ## Find L1R pan file
    bn = os.path.basename(ncf)
    if '_L1R.nc' in bn:
        ncfp = ncf.replace('_L1R.nc', '_L1R_pan.nc')
    if '_L2R.nc' in bn:
        ncfp = ncf.replace('_L2R.nc', '_L1R_pan.nc')

    if 'SEVIRI' in bn:
        ncfp = ncfp.replace('_pan', '_HRV')

    if not os.path.exists(ncfp):
        print('No L1R_pan.nc file available for {}'.format(ncf))
        return

    ## Open gem
    gem = ac.gem.gem(ncf)

    ## parse settings
    if 'user' not in ac.settings:
        ac.settings['user'] = ac.acolite.settings.parse(None, settings=settings, merge=False)
        for k in ac.settings['user']: ac.settings['run'][k] = ac.settings['user'][k]
    setu = ac.acolite.settings.parse(gem.gatts['sensor'], settings=settings)
    for k in ac.settings['user']: setu[k] = ac.settings['user'][k]

    if setu['pans_method'] not in ['panr', 'visr']:
        print('pans_method={} not configured, using panr')
        setu['pans_method'] = 'panr'

    if gem.gatts['sensor'] not in setu['pans_sensors']:
        print('Pan sharpening not supported for {}'.format(gem.gatts['sensor']))
        gem.close()
        return

    ## output file
    ncfo = ncf.replace('.nc', '_pans_{}.nc'.format(setu['pans_method']))
    if output is not None: ncfo = '{}/{}'.format(output, os.path.basename(ncfo))
    print('Pan sharpening {} to {}'.format(ncf, ncfo))
    print('Pan file {}'.format(ncfp))

    # read pan gem
    gemp = ac.gem.gem(ncfp)
    dsp = [ds for ds in gemp.datasets if ('rhot_' in ds) | ('pan' == ds) | ('hrv' == ds)][0]
    rhos_datasets = [ds for ds in gem.datasets if 'rhos_' in ds]
    rhos_wave = np.asarray([float(ds.split('_')[1]) for ds in rhos_datasets])
    rhot_datasets = [ds for ds in gem.datasets if 'rhot_' in ds]
    rhot_wave = np.asarray([float(ds.split('_')[1]) for ds in rhot_datasets])

    ## find rgb datasets
    bgr_ds_rhot = []
    bgr_ds_rhos = []
    for i, wv in enumerate(setu['pans_bgr']):
        wi = np.argsort(np.abs(rhot_wave-float(wv)))[0]
        bgr_ds_rhot.append('rhot_{:.0f}'.format(rhot_wave[wi]))
        bgr_ds_rhos.append('rhos_{:.0f}'.format(rhot_wave[wi]))

    pans_ds = bgr_ds_rhot + bgr_ds_rhos
    if setu['pans_output'] == 'all': pans_ds = [ds for ds in gem.datasets if ds[0:3] == 'rho']

    # get projection info
    nc_projection_ms = gem.nc_projection
    nc_projection_pan = gemp.nc_projection
    if (nc_projection_ms is not None) & (nc_projection_pan is not None):
        fac_x = np.diff(nc_projection_ms['x']['data'])[0] / np.diff(nc_projection_pan['x']['data'])[0]
        fac_y = np.diff(nc_projection_ms['y']['data'])[0] / np.diff(nc_projection_pan['y']['data'])[0]
        factor = fac_x
        print('Assuming pan scale factor {} (x={}, y={})'.format(factor, fac_x, fac_y))
    else:
        factor = setu['pans_scale_factor']
        print('Assuming pan scale factor {}'.format(factor))

    ## read pan data
    pan = gemp.data(dsp)
    gemp.close()
    print(pan.shape)

    print('Using pans_method={}'.format(setu['pans_method']))
    ## compute pan factor
    if setu['pans_method'] == 'panr':
        if dsp in gem.datasets:
            pan_ms = gem.data(dsp)
        else:
            pan_ms = scipy.ndimage.zoom(pan, zoom=1/factor, order=1)

        print(pan_ms.shape)

        pan_ms_ = scipy.ndimage.zoom(pan_ms, factor, order=0)
        print(pan_ms_.shape)

        pan_i = pan_ms_ / pan
    elif setu['pans_method'] == 'visr':
        for i, ds in enumerate(bgr_ds_rhot):
            print(ds)
            bd = gem.data(dsp)
            if i == 0:
                vis_i = bd
            else:
                vis_i += bd
            bd = None
        vis_i = scipy.ndimage.zoom((vis_i)/3, factor, order=0)
        pan_i = vis_i / pan
        vis_i = None
        pan = None

    ## make output file
    gemo = ac.gem.gem(ncfo, new = True)
    gemo.gatts = {g:gem.gatts[g] for g in gem.gatts}
    gemo.nc_projection = nc_projection_pan
    gemo.gatts['ofile'] = ncfo

    ## run through pansharpening
    for ds in gem.datasets:
        if ds in ['transverse_mercator', 'x', 'y']: continue
        if ds not in pans_ds: continue
        data, ds_att = gem.data(ds, attributes = True)
        data_pan = scipy.ndimage.zoom(data, factor, order=0)
        data = None

        if ('rhos_' in ds) | ('rhot_' in ds):
            data_pan /= pan_i
            print('Pan sharpened {}'.format(ds))

        gemo.write(ds, data_pan, ds_att = ds_att)
        print('Wrote {}'.format(ds))
        data_pan = None
    print('Wrote {}'.format(ncfo))
    gem.close()
    gemo.close()

    return(ncfo)
