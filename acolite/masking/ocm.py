## def ocm
## runs OmniCloudMask on ACOLITE L1R/L2R NetCDFs
## optional RGB and RGB-masked and NetCDF outputs
## written by Quinten Vanhellemont, RBINS
## 2025-11-19
## modifications: 2025-12-16 (QV) added dataset keyword, moved wavelengths and rgb range to keywords
##                2026-01-15 (QV) added export_netcdf option, added to acolite.masking
##                2026-02-11 (QV) added masks as keyword

def ocm(ncf, output = None, plot_rgb = True, plot_mask = False, plot_rgb_mask = True,
                  rgb_range = [0, 0.15], rgb_default = [665, 560, 490], wave_ref = [665, 560, 885],
                  fill_nan = True, fill_value = 0, dataset = 'rhot', masks = None,
                  export_maps = False, export_netcdf = False, return_mask = True):
    import os
    import acolite as ac
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    import numpy as np
    from omnicloudmask import predict_from_array

    bn = os.path.basename(ncf)
    odir = '{}'.format(output) if output is not None else os.path.dirname(ncf)
    if not os.path.exists(odir): os.makedirs(odir)

    bno, ext = os.path.splitext(bn)
    gatts = ac.shared.nc_gatts(ncf)
    if 'sensor_lut' in gatts:
        sensor = gatts['sensor_lut']
    else:
        sensor = gatts['sensor']
        setu = ac.acolite.settings.merge(sensor = sensor, settings = None)
        if setu['rsr_version'] is not None:
            sensor += '_{}'.format(setu['rsr_version'])

    rsrd = ac.shared.rsr_dict(sensor = sensor)
    wave_sen = np.asarray([rsrd[sensor]['wave_nm'][b] for b in rsrd[sensor]['rsr_bands']])

    ## cloud mask rgb, green, nir
    print('Computing cloud mask')
    wave_idx = [np.argsort(np.abs(wave_sen - w))[0] for w in wave_ref]
    wave_sel = ['{:.0f}'.format(wave_sen[i]) for i in wave_idx]

    for wi, w in enumerate(wave_sel):
        ds = '{}_{}'.format(dataset, wave_sel[wi])
        print('Reading {}'.format(ds))
        d = ac.shared.nc_data(ncf, ds)
        if fill_nan:
            d[np.isnan(d)] = fill_value
            d[d.mask] = fill_value

        if wi == 0:
            rhot = d
        else:
            rhot = np.dstack((rhot, d))
    ## rearrange axes for OmniCloudMask
    input_array = np.moveaxis(rhot, (2), (0))
    pred_mask = predict_from_array(input_array)[0, :, :]

    ## export NetCDF file
    if export_netcdf:
        if output is None:
            odir = os.path.dirname(os.path.abspath(ncf))
        else:
            odir = os.path.abspath(output)
        bn = os.path.basename(ncf)
        ofile = '{}/{}_OCM.nc'.format(odir, bn[0:-3])
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {g: gatts[g] for g in gatts}
        gemo.nc_projection = ac.shared.nc_projection_read(ncf)
        gemo.write('lat', ac.shared.nc_data(ncf, 'lat'))
        gemo.write('lon', ac.shared.nc_data(ncf, 'lon'))
        gemo.write('ocm_pred_mask', pred_mask)
        gemo.close()
        print('Wrote {}'.format(ofile))
    ## end export NetCDF file

    ## create rgb/maps
    if export_maps:
        print('Creating RGB')
        wave_idx = [np.argsort(np.abs(wave_sen - w))[0] for w in rgb_default]
        wave_rgb = ['{:.0f}'.format(wave_sen[i]) for i in wave_idx]
        for wi, w in enumerate(wave_rgb):
            ds = '{}_{}'.format(dataset, wave_rgb[wi])
            print('Reading {}'.format(ds))
            d = ac.shared.nc_data(ncf, ds)
            if wi == 0:
                rgb = d
            else:
                rgb = np.dstack((rgb, d))

        rgb_scaled = ac.shared.datascl(rgb, dmin = rgb_range[0], dmax = rgb_range[1])

        if plot_rgb:
            rgb_out = '{}/{}_{}_RGB.png'.format(odir, bno, dataset)
            plt.imshow(rgb_scaled)
            plt.title('{} {} \n RGB'.format(sensor.replace('_', '/'), gatts['isodate'][0:19]))
            plt.axis('off')
            plt.savefig(rgb_out, dpi = 300, bbox_inches = 'tight', facecolor = 'white')
            plt.close()
            print('Wrote RGB to {}'.format(rgb_out))

        ## colour masks
        if masks is None:
            masks = {1: {'name': 'Thick Cloud', 'color': 'red'},
                     2: {'name': 'Thin Cloud', 'color': 'orange'},
                     3: {'name': 'Cloud Shadow', 'color': 'yellow'}}

        ## make mask figure
        if plot_mask:
            mask_out = '{}/{}_{}_mask.png'.format(odir, bno, dataset)
            for m in masks:
                cmap = mpl.colors.ListedColormap(mpl.colors.to_rgb(masks[m]['color']))
                mask = (pred_mask == m) * 1.0
                mask[mask == 0] = np.nan
                plt.imshow(mask, cmap = cmap)
            plt.title('{} {} \n OmniCloudMask'.format(sensor.replace('_', '/'), gatts['isodate'][0:19]))
            plt.axis('off')
            plt.savefig(mask_out, dpi = 300, bbox_inches = 'tight', facecolor = 'white')
            plt.close()
            print('Wrote cloud mask to {}'.format(mask_out))

        ## make RGB mask figure
        if plot_rgb_mask:
            rgb_mask_out = '{}/{}_{}_RGB_mask.png'.format(odir, bno, dataset)

            plt.imshow(rgb_scaled)
            xlim = plt.xlim()
            ylim = plt.ylim()

            for m in masks:
                cmap = mpl.colors.ListedColormap(mpl.colors.to_rgb(masks[m]['color']))
                mask = (pred_mask == m) * 1.0
                mask[mask == 0] = np.nan
                plt.imshow(mask, cmap = cmap)

                r = plt.scatter(-10,-10, color = masks[m]['color'], linestyle = '', label = masks[m]['name'], marker = 's')

            plt.xlim(xlim)
            plt.ylim(ylim)
            plt.legend()
            plt.title('{} {} \n RGB+OmniCloudMask'.format(sensor.replace('_', '/'), gatts['isodate'][0:19]))
            plt.axis('off')
            plt.savefig(rgb_mask_out, dpi = 300, bbox_inches = 'tight', facecolor = 'white')
            plt.close()
            print('Wrote cloud masked RGB to {}'.format(rgb_mask_out))
    ## end create rgb/maps

    if return_mask:
        return(pred_mask)
