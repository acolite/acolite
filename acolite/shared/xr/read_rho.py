## function to read reflectance data from NetCDF using xarray
##
## written by Quinten Vanhellemont, RBINS
## 2026-06-10
## modifications:

def read_rho(ncf, parameter = 'rhot', datasets = None):
    import xarray as xr
    import numpy as np

    ## additional data to read
    if datasets is not None:
        if type(datasets) is not list:
            datasets = [datasets]

    ## set up empty lists
    wavelength = []
    width = []
    rho = []
    data = {}

    ## open NetCDF
    with xr.open_dataset(ncf) as ds:
        nc_datasets = list(ds.keys())

        ## read attributes
        gatts = {k: ds.attrs[k] for k in ds.attrs}

        ## read rhot
        for d in nc_datasets:
            if d.startswith('{}_'.format(parameter)):
                wv = None if 'wave_nm' not in ds[d].attrs else ds[d].attrs['wave_nm']
                wavelength.append(wv)
                wd = None if 'width' not in ds[d].attrs else ds[d].attrs['width']
                width.append(wd)
                rho.append(ds[d].values)

        ## read additional dataset
        if type(datasets) == list:
            for d in datasets:
                if d in ds:
                    data[d] = ds[d].values

    ## convert to arrays
    data[parameter] = np.asarray(rho)
    del rho
    data['wavelength'] = np.asarray(wavelength)
    del wavelength
    data['width'] = np.asarray(width)
    del width

    return(data, gatts)
