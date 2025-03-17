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

## read global attributes, and all datasets from NetCDF (with/without group)
def nc_read_all(ncf, group = None):
    from netCDF4 import Dataset
    with Dataset(ncf) as nc:
        gatts = {attr : getattr(nc,attr) for attr in nc.ncattrs()}
        if group is not None:
            if group in nc.groups: nc = nc.groups[group]
        datasets = nc.variables.keys()
        data = {}
        atts = {}
        for ds in datasets:
            data[ds] = nc.variables[ds][:]
            atts[ds] = {a : getattr(nc.variables[ds],a) for a in nc.variables[ds].ncattrs()}
    return(gatts, data, atts)

# read dataset from netcdf
# Last updates: 2016-12-19 (QV) added crop (x0,x1,y0,y1)
##              2017-03-16 (QV) added sub keyword (xoff, yoff, xcount, ycount)
##              2023-10-31 (QV) added dtype keyword
##              2024-02-21 (QV) assume 3d array subset is keeps the first dimension
##              2024-07-01 (QV) add index for 3d array subset

def nc_data(file, dataset, crop=False, sub=None, attributes=False, group=None, dtype=None, axis_3d = 0):
    from netCDF4 import Dataset
    with Dataset(file) as nc:
        if group is not None:
            if group in nc.groups: nc = nc.groups[group]
        if sub is None:
            if crop is False:
                data = nc.variables[dataset][:]
            else:
                if len(crop) == 4:
                    if len(nc.variables[dataset].shape) == 3:
                        if axis_3d == 0:
                            data = nc.variables[dataset][:, crop[2]:crop[3]:1,crop[0]:crop[1]:1]
                        if axis_3d == 1:
                            data = nc.variables[dataset][crop[2]:crop[3]:1, :, crop[0]:crop[1]:1]
                        if axis_3d == 2:
                            data = nc.variables[dataset][crop[2]:crop[3]:1,crop[0]:crop[1]:1, :]
                    else:
                        data = nc.variables[dataset][crop[2]:crop[3]:1,crop[0]:crop[1]:1]
                else: data = nc.variables[dataset][:]
        else:
            if len(sub) == 4:
                if len(nc.variables[dataset].shape) == 3:
                    if axis_3d == 0:
                        data = nc.variables[dataset][:, sub[1]:sub[1]+sub[3]:1,sub[0]:sub[0]+sub[2]:1]
                    if axis_3d == 1:
                        data = nc.variables[dataset][sub[1]:sub[1]+sub[3]:1,:,sub[0]:sub[0]+sub[2]:1]
                    if axis_3d == 2:
                        data = nc.variables[dataset][sub[1]:sub[1]+sub[3]:1,sub[0]:sub[0]+sub[2]:1,:]
                else:
                    data = nc.variables[dataset][sub[1]:sub[1]+sub[3]:1,sub[0]:sub[0]+sub[2]:1]
            else: data = nc.variables[dataset][:]
        if attributes:
            atts = {attr : getattr(nc.variables[dataset],attr) for attr in nc.variables[dataset].ncattrs()}
    if dtype is not None: data = data.astype(dtype)

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

# read groups from netcdf
def nc_groups(file, datasets = True):
    from netCDF4 import Dataset
    with Dataset(file) as nc:
        if datasets:
            groups = {}
            for g in nc.groups:
                groups[g] = [v for v in nc.groups[g].variables]
        else:
            groups = [g for g in nc.groups]
    return groups

# read datasets in netcdf
# Last updates: 2023-05-31 (QV) added group keyword
def nc_datasets(file, group = None):
    from netCDF4 import Dataset
    with Dataset(file) as nc:
        if group is not None:
            if group in nc.groups: nc = nc.groups[group]
        ds = list(nc.variables.keys())
    return ds
