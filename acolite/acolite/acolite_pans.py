## def acolite_pans
## function to pansharpen ACOLITE NetCDF outputs
## written by Quinten Vanhellemont, RBINS
## 2022-02-20
## modifications: 2022-02-21 (QV) acolite integration

def acolite_pans(ncf, output = None, settings = {}):

    import os
    import acolite as ac
    import numpy as np
    import scipy.ndimage

    ## Find L1R pan file
    ncfp = ncf.replace('_L2R.nc', '_L1R_pan.nc')
    if not os.path.exists(ncfp):
        print('No L1R_pan.nc file available for {}'.format(ncf))
        return()

    ## parse settings
    setu = ac.acolite.settings.parse(None, settings=settings)
    if setu['pans_method'] not in ['panr', 'visr']:
        print('pans_method={} not configured, using panr')
        setu['pans_method'] = 'panr'

    ## Open gatts
    gatts = ac.shared.nc_gatts(ncf)
    if gatts['sensor'] not in setu['pans_sensors']:
        print('Pan sharpening not supported for {}'.format(gatts['sensor']))
        return()

    ## output file
    ncfo = ncf.replace('.nc', '_pans_{}.nc'.format(setu['pans_method']))
    if output is not None: ncfo = '{}/{}'.format(output, os.path.basename(ncfo))
    print('Pan sharpening {} to {}'.format(ncf, ncfo))

    # read pan datasets
    datasets_pan = ac.shared.nc_datasets(ncfp)
    dsp = [ds for ds in datasets_pan if 'rhot_' in ds][0]

    # read ms datasets
    datasets = ac.shared.nc_datasets(ncf)
    rhos_datasets = [ds for ds in datasets if 'rhos_' in ds]
    rhos_wave = np.asarray([float(ds.split('_')[1]) for ds in rhos_datasets])
    rhot_datasets = [ds for ds in datasets if 'rhot_' in ds]
    rhot_wave = np.asarray([float(ds.split('_')[1]) for ds in rhot_datasets])

    ## find rgb datasets
    bgr_ds_rhot = []
    bgr_ds_rhos = []
    for i, wv in enumerate(setu['pans_bgr']):
        wi = np.argsort(np.abs(rhot_wave-float(wv)))[0]
        bgr_ds_rhot.append('rhot_{:.0f}'.format(rhot_wave[wi]))
        bgr_ds_rhos.append('rhos_{:.0f}'.format(rhot_wave[wi]))

    pans_ds = bgr_ds_rhot + bgr_ds_rhos
    if setu['pans_output'] == 'all': pans_ds = [ds for ds in datasets if ds[0:3] == 'rho']

    # get projection info
    nc_projection_ms = ac.shared.nc_read_projection(ncf)
    nc_projection_pan = ac.shared.nc_read_projection(ncfp)
    fac_x = np.diff(nc_projection_ms['x']['data'])[0] / np.diff(nc_projection_pan['x']['data'])[0]
    fac_y = np.diff(nc_projection_ms['y']['data'])[0] / np.diff(nc_projection_pan['y']['data'])[0]
    factor = fac_x
    print('Assuming pan scale factor {} (x={}, y={})'.format(factor, fac_x, fac_y))

    ## read pan data
    pan = ac.shared.nc_data(ncfp, dsp)

    print('Using pans_method={}'.format(setu['pans_method']))
    ## compute pan factor
    if setu['pans_method'] == 'panr':
        pan_ms = ac.shared.nc_data(ncf, dsp)
        pan_ms_ = scipy.ndimage.zoom(pan_ms, factor, order=0)
        pan_i = pan_ms_ / pan
    elif setu['pans_method'] == 'visr':
        for i, ds in enumerate(bgr_ds_rhot):
            print(ds)
            bd = ac.shared.nc_data(ncf, ds)
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
    gatts_out = {g:gatts[g] for g in gatts}
    gatts_out['ofile'] = ncfo

    ## run through pansharpening
    new = True
    for ds in datasets:
        if ds in ['transverse_mercator', 'x', 'y']: continue
        if ds not in pans_ds: continue

        data, ds_att = ac.shared.nc_data(ncf, ds, attributes=True)
        data_pan = scipy.ndimage.zoom(data, factor, order=0)
        data = None

        if ('rhos_' in ds) | ('rhot_' in ds):
            data_pan /= pan_i
            print('Pan sharpened {}'.format(ds))

        ac.output.nc_write(ncfo, ds, data_pan, dataset_attributes=ds_att,
                           attributes=gatts_out, nc_projection=nc_projection_pan, new=new)
        print('Wrote {}'.format(ds))
        data_pan = None
        new = False
    print('Wrote {}'.format(ncfo))

    ## mapping - could be moved to acolite_run
    setu_map = {s:setu[s] for s in setu}
    setu_map['rgb_rhot'] = setu_map['pans_rgb_rhot']
    setu_map['rgb_rhos'] = setu_map['pans_rgb_rhos']
    if setu_map['pans_rgb_rhot'] or setu_map['pans_rgb_rhos']:
        ac.acolite.acolite_map(ncfo, plot_all=False, settings=setu_map)
    if setu['pans_export_geotiff_rgb']:
        ac.output.nc_to_geotiff_rgb(ncfo, settings=setu_map)

    return(ncfo)
