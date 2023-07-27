## series of NetCDF reading functions
## written by Quinten Vanhellemont, RBINS
##
## modifications:
##              2021-11-17 (QV) updated file handling
##              2023-05-31 (QV) added group keyword

# read dataset and global attributes from netcdf
def nc_read(file, dataset, group = None):
    from netCDF4 import Dataset
    with Dataset(file) as nc:
        if group is not None:
            if group in nc.groups: nc = nc.groups[group]
        gatts = {attr : getattr(nc,attr) for attr in nc.ncattrs()}
        out_array = nc.variables[dataset][:]
    return (out_array, gatts)

# read dataset from netcdf
# Last updates: 2016-12-19 (QV) added crop (x0,x1,y0,y1)
##              2017-03-16 (QV) added sub keyword (xoff, yoff, xcount, ycount)
def nc_data(file, dataset, crop=False, sub=None, attributes=False, group=None):
    from netCDF4 import Dataset
    with Dataset(file) as nc:
        if group is not None:
            if group in nc.groups: nc = nc.groups[group]
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
def nc_atts(file, dataset, group = None):
    from netCDF4 import Dataset
    with Dataset(file) as nc:
        if group is not None:
            if group in nc.groups: nc = nc.groups[group]
        atts = {attr : getattr(nc.variables[dataset],attr) for attr in nc.variables[dataset].ncattrs()}
    return atts

# read global attributes from netcdf
def nc_gatts(file):
    from netCDF4 import Dataset
    with Dataset(file) as nc:
        gatts = {attr : getattr(nc,attr) for attr in nc.ncattrs()}
    return gatts

# read datasets in netcdf
# Last updates: 2023-05-31 (QV) added group keyword
def nc_datasets(file, group = None):
    from netCDF4 import Dataset
    with Dataset(file) as nc:
        if group is not None:
            if group in nc.groups: nc = nc.groups[group]
        ds = list(nc.variables.keys())
    return ds
