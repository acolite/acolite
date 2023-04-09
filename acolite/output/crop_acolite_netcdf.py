## def crop_acolite_netcdf
## sample script for cropping acolite NetCDF file
## written by Quinten Vanhellemont, RBINS
## 2023-04-02
## modifications: 2023-04-09 (QV) Fix for NetCDF without nc_projection

def crop_acolite_netcdf(ncf, output=None, limit=None, polygon=None):
    import os
    import acolite as ac

    if polygon is not None:
        limit = ac.shared.polygon_limit(polygon)

    if limit is None:
        print("Please provide a four element limit [S, W, N, E] or a polygon file for L1R cropping")
        return

    gatts = ac.shared.nc_gatts(ncf)
    dn = os.path.dirname(ncf)
    bn = os.path.basename(ncf)
    if output is None: output = dn

    datasets = ac.shared.nc_datasets(ncf)
    nc_projection = ac.shared.nc_read_projection(ncf)

    ofile = '{}/{}'.format(dn, bn.replace('.nc', '_crop.nc'))

    ## load lat/lon and crop
    lat, latatt = ac.shared.nc_data(ncf, 'lat', attributes=True)
    lon, lonatt = ac.shared.nc_data(ncf, 'lon', attributes=True)
    sub = ac.shared.geolocation_sub(lat, lon, limit)

    if sub is not None:
        print('Cropping file {} to {}'.format(ncf,ofile))
        print('New limit {}'.format(limit))
        print('Sub {}'.format(sub))
        if nc_projection is not None:
            nc_projection['x']['data'] = nc_projection['x']['data'][sub[0]:sub[0]+sub[2]]
            nc_projection['y']['data'] = nc_projection['y']['data'][sub[1]:sub[1]+sub[3]]

        new = True
        for ds in datasets:
            if nc_projection is not None:
                if ds in nc_projection: continue
            d, da = ac.shared.nc_data(ncf, ds, attributes=True, sub=sub)

            ## output dataset
            ac.output.nc_write(ofile, ds, d, dataset_attributes = da,
                               new=new, attributes=gatts, nc_projection=nc_projection)
            print('Wrote {} ({}) to {}'.format(ds, d.shape, ofile))
            new = False
        return(ofile)
