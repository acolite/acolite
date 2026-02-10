## def crop_acolite_netcdf
## sample script for cropping acolite NetCDF file
## written by Quinten Vanhellemont, RBINS
## 2023-04-02
## modifications: 2023-04-09 (QV) Fix for NetCDF without nc_projection
##                2024-04-16 (QV) Use gem for input/output
##                2025-01-30 (QV) removed polygon keyword
##                2026-02-10 (QV) added sub, crop, and crop_name keywords

def crop_acolite_netcdf(ncf, output = None, limit = None,
                        sub = None, crop = None, crop_name = 'crop',
                        verbosity = 0):
    import os
    import acolite as ac

    if (limit is None) & (sub is None) & (crop is None):
        print("Please provide a four element limit [S, W, N, E] or subset [x0, nx, y0, ny] or crop subset [x0, x1, y0, y1] for L1R cropping")
        return

    gem = ac.gem.gem(ncf)
    dn = os.path.dirname(ncf)
    bn = os.path.basename(ncf)
    if output is None: output = dn

    ofile = '{}/{}'.format(output, bn.replace('.nc', '_{}.nc'.format(crop_name)))

    ## load lat/lon and crop
    if (sub is None) & (crop is None):
        lat, latatt = gem.data('lat', attributes = True)
        lon, lonatt = gem.data('lon', attributes = True)
        sub = ac.shared.geolocation_sub(lat, lon, limit)

    if (sub is None):
        sub = crop[0], crop[2], crop[1]-crop[0], crop[3]-crop[2]

    if sub is not None:
        if verbosity > 1:
            print('Cropping file {} to {}'.format(ncf,ofile))
            print('New limit {}'.format(limit))
            print('Sub {}'.format(sub))

        ## setup output dataset
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gem.gatts[k] for k in gem.gatts} ## set output attributes
        gemo.gatts['ofile'] = ofile
        gemo.verbosity = verbosity

        if gem.nc_projection is not None:
            gemo.nc_projection = gem.nc_projection
            gemo.nc_projection['x']['data'] = gemo.nc_projection['x']['data'][sub[0]:sub[0]+sub[2]]
            gemo.nc_projection['y']['data'] = gemo.nc_projection['y']['data'][sub[1]:sub[1]+sub[3]]

        new = True
        for ds in gem.datasets:
            if gem.nc_projection is not None:
                if ds in gem.nc_projection_keys: continue

            try:
                d, da = gem.data(ds, attributes = True, sub = sub)
            except:
                if verbosity > 0:
                    print('Could not crop dataset {}'.format(ds))
                continue

            ## output dataset
            gemo.write(ds, d, ds_att = da)
            if verbosity > 1:
                print('Wrote {} ({}) to {}'.format(ds, d.shape, ofile))
        return(ofile)
