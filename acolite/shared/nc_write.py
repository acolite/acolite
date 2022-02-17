## QV 2021-11-27 added nc_gatts_update

# update global attributes from netcdf
def nc_gatts_update(file, gatts):
    from netCDF4 import Dataset
    with Dataset(file, 'a', format='NETCDF4') as nc:
        for key in gatts.keys():
            if gatts[key] is not None:
                try:
                    setattr(nc, key, gatts[key])
                except:
                    print('Failed to write attribute: {}'.format(key))
