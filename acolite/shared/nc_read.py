## QV 2021-11-17 updated file handling

# read dataset and global attributes from netcdf
def nc_read(file, dataset):
    from netCDF4 import Dataset
    with Dataset(file) as nc:
        gatts = {attr : getattr(nc,attr) for attr in nc.ncattrs()}
        out_array = nc.variables[dataset][:]
    return (out_array, gatts)

# read dataset from netcdf
# Last updates: 2016-12-19 (QV) added crop (x0,x1,y0,y1)
##              2017-03-16 (QV) added sub keyword (xoff, yoff, xcount, ycount)
def nc_data(file, dataset, crop=False, sub=None, attributes=False):
    from netCDF4 import Dataset
    with Dataset(file) as nc:
        if sub is None:
            if crop is False:
                data = nc.variables[dataset][:]
            else:
                if len(crop) == 4: data = nc.variables[dataset][crop[2]:crop[3]:1,crop[0]:crop[1]:1]
                else: data = nc.variables[dataset][:]
        else:
            if len(sub) == 4: data = nc.variables[dataset][sub[1]:sub[1]+sub[3]:1,sub[0]:sub[0]+sub[2]:1]
            else: data = nc.variables[dataset][:]
        if attributes:
            atts = {attr : getattr(nc.variables[dataset],attr) for attr in nc.variables[dataset].ncattrs()}
    if attributes:
        return(data,atts)
    else:
        return(data)

## get attributes for given dataset
def nc_atts(file, dataset):
    from netCDF4 import Dataset
    with Dataset(file) as nc:
        atts = {attr : getattr(nc.variables[dataset],attr) for attr in nc.variables[dataset].ncattrs()}
    return atts

# read global attributes from netcdf
def nc_gatts(file):
    from netCDF4 import Dataset
    with Dataset(file) as nc:
        gatts = {attr : getattr(nc,attr) for attr in nc.ncattrs()}
    return gatts

# read datasets in netcdf
def nc_datasets(file):
    from netCDF4 import Dataset
    with Dataset(file) as nc:
        ds = list(nc.variables.keys())
    return ds
