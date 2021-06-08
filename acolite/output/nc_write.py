## def write_rgb
## writes dataset to netcdf file
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07
## modifications: QV 2017-01-18 converted to function
##                QV 2017-01-25 added attributes keyword to write global attributes to new file
##                QV 2017-04-18 convert wavelength to float to avoid error in BEAM/SNAP
##                QV 2017-04-20 added data type check for ndarrays, check for datasets already in file
##                QV 2017-06-21 added chunking option/setting
##                QV 2017-11-27 added offset option, and replace_nan option
##                QV 2017-12-05 made "metadata" a keyword; added dataset_attributes keyword
##                QV 2017-12-06 added global_dims keyword
##                QV 2018-03-14 added nan writing for new datasets with offset
##                QV 2018-04-17 changed to float32 for float datasets, added Rrs to auto_grouping
##                QV 2018-07-18 changed datatype for writing, to avoid int overflow
##                QV 2018-07-24 changed global attributes
##                QV 2020-07-14 added fillvalue keyword
##                QV 2020-07-22 added update_attributes keyword
##                QV 2020-07-23 skip fillvalue in ds attributes
##                QV 2021-02-09 added replace_nan option for data without offset, changed numpy import
##                QV 2021-06-04 added dataset attributes defaults

def nc_write(ncfile, dataset, data, wavelength=None, global_dims=None,
                 new=False, attributes=None, update_attributes=False,
                 keep=True, offset=None, replace_nan=False, metadata=None, dataset_attributes=None, double=False,
                 chunking=True, chunk_tiles=[10,10], chunksizes=None, fillvalue=None,
                 format='NETCDF4',#'NETCDF4_CLASSIC',
                 nc_compression=False # currently off: file saving takes *much* longer,
                                      #about 30% file size reduction for Pléiades
                 ):

    from netCDF4 import Dataset
    import time, os
    from math import ceil
    import numpy as np

    import re
    import acolite as ac

    ## import atts for dataset
    atts = None
    for p in ac.param['attributes']:
        if p['parameter'] == dataset: atts = {t:p[t] for t in p}
    if atts is None:
        for p in ac.param['attributes']:
            if re.match(p['parameter'], dataset):
                atts = {t:p[t] for t in p}
                try:
                    wave = int(dataset.split('_')[-1])
                    atts['wavelength'] = wave
                except:
                    pass

    ## set attributes from provided/defaults
    if atts is not None:
        if dataset_attributes is None:
            dataset_attributes = {t:atts[t] for t in atts}
        else:
            for t in atts:
                if t not in dataset_attributes:
                    dataset_attributes[t] = atts[t]
        dataset_attributes['parameter'] = dataset

    if os.path.exists(os.path.dirname(ncfile)) is False:
         os.makedirs(os.path.dirname(ncfile))

    dims = data.shape
    if global_dims is None: global_dims = dims

    if chunking:
        if chunksizes is not None:
            chunksizes=(ceil(dim[0]/chunk_tiles[0]), ceil(dim[1]/chunk_tiles[1]))

    if new:
        if os.path.exists(ncfile): os.remove(ncfile)
        #nc = Dataset(ncfile, 'w', format='NETCDF4_CLASSIC')
        nc = Dataset(ncfile, 'w', format=format)

        ## set global attributes
        setattr(nc, 'generated_by', 'ACOLITE' )
        setattr(nc, 'generated_on',time.strftime('%Y-%m-%d %H:%M:%S %Z'))
        #setattr(nc, 'project', 'PONDER' )
        setattr(nc, 'contact', 'Quinten Vanhellemont' )

        ## set beam dataformat global attributes
        setattr(nc, 'product_type', 'NetCDF' )
        setattr(nc, 'metadata_profile', 'beam' )
        setattr(nc, 'metadata_version', '0.5' )
        #setattr(nc, 'auto_grouping', 'rtoa:rrc:rhos' )
        setattr(nc, 'auto_grouping', 'rhot:rhorc:rhos:rhow:Rrs')

        ## CF convention
        setattr(nc, 'Conventions', 'CF-1.7')
        ## to add: title , history , institution , source , comment and references

        if attributes is not None:
            for key in attributes.keys():
                if attributes[key] is not None:
                    try:
                        setattr(nc, key, attributes[key])
                    except:
                        print('Failed to write attribute: {}'.format(key))

        ## set up x and y dimensions
        nc.createDimension('x', global_dims[1])
        nc.createDimension('y', global_dims[0])
    else:
        nc = Dataset(ncfile, 'a', format=format)
        if update_attributes:
            if attributes is not None:
                for key in attributes.keys():
                    if attributes[key] is not None:
                        try:
                            setattr(nc, key, attributes[key])
                        except:
                            print('Failed to write attribute: {}'.format(key))

    if (not double) & (data.dtype == np.float64):
        data = data.astype(np.float32)

    ## write data
    if dataset in nc.variables.keys():
        ## dataset already in NC file
        if offset is None:
            if replace_nan:
                tmp = nc.variables[dataset][:]
                sub_isnan=np.where(np.isnan(tmp))
                tmp[sub_isnan]=data[sub_isnan]
                nc.variables[dataset][:] = tmp
                #nc.variables[dataset][sub_isnan] = data[sub_isnan]
                tmp = None
            else:
                if data.dtype in (np.float32, np.float64): nc.variables[dataset][:] = np.nan
                nc.variables[dataset][:] = data
        else:
            if replace_nan:
                tmp = nc.variables[dataset][offset[1]:offset[1]+dims[0],offset[0]:offset[0]+dims[1]]
                sub_isnan=np.where(np.isnan(tmp))
                tmp[sub_isnan]=data[sub_isnan]
                nc.variables[dataset][offset[1]:offset[1]+dims[0],offset[0]:offset[0]+dims[1]] = tmp
                tmp = None
            else:
                nc.variables[dataset][offset[1]:offset[1]+dims[0],offset[0]:offset[0]+dims[1]] = data
    else:
        ## new dataset
        var = nc.createVariable(dataset,data.dtype,('y','x'),
                                fill_value=fillvalue,
                                zlib=nc_compression,chunksizes=chunksizes)
        if wavelength is not None: setattr(var, 'wavelength', float(wavelength))
        ## set attributes
        if dataset_attributes is not None:
            for att in dataset_attributes.keys():
                if att in ['_FillValue']: continue
                setattr(var, att, dataset_attributes[att])

        if offset is None:
            if data.dtype in (np.float32, np.float64): var[:] = np.nan
            var[:] = data
        else:
            if data.dtype in (np.float32, np.float64): var[:] = np.nan
            var[offset[1]:offset[1]+dims[0],offset[0]:offset[0]+dims[1]] = data
    if keep is not True: data = None

    ## close netcdf file
    nc.close()
    nc=None
