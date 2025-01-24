## def nc_write
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
##                QV 2021-07-19 change to using setncattr
##                QV 2021-12-08 added nc_projection
##                QV 2022-11-09 added option to update nc_projection
##                QV 2023-07-12 added discretisation
##                QV 2024-01-31 added skip_attributes
##                2024-04-15 (QV) allow passing of netCDF4 Dataset
##                2024-04-16 (QV) removed NetCDF compression parameters from keywords
##                2025-01-23 (QV) added break to dataset attributes check

def nc_write(ncfile, dataset, data, wavelength=None, global_dims=None,
                 new=False, attributes=None, update_attributes=False,
                 keep=True, offset=None, replace_nan=False,
                 metadata=None, dataset_attributes=None, double=False,
                 chunking=True, chunk_tiles=[10,10], chunksizes=None, fillvalue=None,
                 nc_projection = None, update_projection=False,
                 format='NETCDF4', return_nc = False):


    from netCDF4 import Dataset
    import time, os
    from math import ceil
    import numpy as np

    import re
    import acolite as ac

    open_file = True
    if type(ncfile) is Dataset:
        open_file, new = False, False

    ## import atts for dataset
    atts = None
    for p in ac.param['attributes']:
        if p['parameter'] == dataset:
            atts = {t:p[t] for t in p}
            break
    if atts is None:
        for p in ac.param['attributes']:
            if re.match(p['parameter'], dataset) is not None:
                atts = {t:p[t] for t in p}
                if p['parameter'][0:2] != 'bt':
                    try:
                        wave = int(dataset.split('_')[-1])
                        atts['wavelength'] = wave
                    except:
                        pass
                break

    ## load discretisation settings
    pdisc, discretise = None, False
    netcdf_discretisation = ac.settings['run']['netcdf_discretisation']
    if netcdf_discretisation:
        if dataset in ac.param['discretisation']:
            pdisc = ac.param['discretisation'][dataset]
        else:
            for p in ac.param['discretisation']:
                if re.match(p, dataset): pdisc = ac.param['discretisation'][p]
    if pdisc is not None: discretise = pdisc['discretise']

    ## set fill value to minimum in the discretisation
    if discretise:
        if pdisc['source_type'] != pdisc['target_type']:
            fillvalue = pdisc['add_offset']

    ## compression options for current run
    netcdf_compression = ac.settings['run']['netcdf_compression']
    netcdf_compression_level = ac.settings['run']['netcdf_compression_level']
    netcdf_compression_least_significant_digit = ac.settings['run']['netcdf_compression_least_significant_digit']

    ## set attributes from provided/defaults
    if atts is not None:
        if dataset_attributes is None:
            dataset_attributes = {t:atts[t] for t in atts}
        else:
            for t in atts:
                if t not in dataset_attributes:
                    dataset_attributes[t] = atts[t]
        dataset_attributes['parameter'] = dataset
    ## convert bool attributes to string
    if dataset_attributes is not None:
        for att in dataset_attributes:
            if type(dataset_attributes[att]) == bool:
                dataset_attributes[att] = str(dataset_attributes[att])

    dims = data.shape
    if global_dims is None: global_dims = dims

    if chunking:
        if chunksizes is not None:
            chunksizes=(ceil(dim[0]/chunk_tiles[0]), ceil(dim[1]/chunk_tiles[1]))

    ## create new file
    if new:
        ## remove existing file
        if os.path.exists(ncfile): os.remove(ncfile)
        ## make output directory
        if not os.path.exists(os.path.dirname(ncfile)): os.makedirs(os.path.dirname(ncfile))

        nc = Dataset(ncfile, 'w', format=format)

        ## set global attributes
        nc.setncattr('generated_by', 'ACOLITE' )
        nc.setncattr('generated_on',time.strftime('%Y-%m-%d %H:%M:%S %Z'))
        nc.setncattr('contact', 'Quinten Vanhellemont' )
        nc.setncattr('acolite_version', ac.version )

        ## set beam dataformat global attributes
        nc.setncattr('product_type', 'NetCDF' )
        nc.setncattr('metadata_profile', 'beam' )
        nc.setncattr('metadata_version', '0.5' )
        nc.setncattr('auto_grouping', 'rhot:rhorc:rhos:rhow:Rrs:Lt:Ed')

        ## CF convention
        nc.setncattr('Conventions', 'CF-1.7')
        ## to add: title , history , institution , source , comment and references

        if attributes is not None:
            for key in attributes.keys():
                if key in ac.config['skip_attributes']: continue
                if attributes[key] is not None:
                    try:
                        nc.setncattr(key, attributes[key])
                    except:
                        print('Failed to write attribute: {}'.format(key))

        ## set up x and y dimensions
        x = nc.createDimension('x', global_dims[1])
        y = nc.createDimension('y', global_dims[0])

        ## set up NetCDF projection if provided
        if nc_projection is not None:
            pkey = [k for k in nc_projection.keys() if k not in ['x', 'y']][0]
            nc.setncattr('projection_key', pkey)

            var = nc.createVariable(pkey, np.float64)
            for att in nc_projection[pkey]['attributes'].keys():
                if att in ['_FillValue']: continue
                var.setncattr(att, nc_projection[pkey]['attributes'][att])

            for v in ['x', 'y']:
                var = nc.createVariable(v, nc_projection[v]['data'].dtype, (v,))
                #var[:] = nc_projection[v]['data']
                var = nc.variables[v]
                var[:] = nc_projection[v]['data']
                for att in nc_projection[v]['attributes'].keys():
                    if att in ['_FillValue']: continue
                    var.setncattr(att, nc_projection[v]['attributes'][att])
    else:
        if open_file:
            nc = Dataset(ncfile, 'a', format=format)
        else:
            nc = ncfile

        if update_attributes:
            if attributes is not None:
                for key in attributes.keys():
                    if key in ac.config['skip_attributes']: continue
                    if attributes[key] is not None:
                        try:
                            nc.setncattr(key, attributes[key])
                        except:
                            print('Failed to write attribute: {}'.format(key))

        ## set up NetCDF projection if provided
        if (update_projection) & (nc_projection is not None):
            pkey = [k for k in nc_projection.keys() if k not in ['x', 'y']][0]
            nc.setncattr('projection_key', pkey)

            var = nc.createVariable(pkey, np.float64)
            for att in nc_projection[pkey]['attributes'].keys():
                if att in ['_FillValue']: continue
                var.setncattr(att, nc_projection[pkey]['attributes'][att])

            for v in ['x', 'y']:
                if v not in nc.variables.keys():
                    var = nc.createVariable(v, nc_projection[v]['data'].dtype, (v,))
                #var[:] = nc_projection[v]['data']
                var = nc.variables[v]
                var[:] = nc_projection[v]['data']
                for att in nc_projection[v]['attributes'].keys():
                    if att in ['_FillValue']: continue
                    var.setncattr(att, nc_projection[v]['attributes'][att])


    if (not double) & (data.dtype == np.float64):
        data = data.astype(np.float32)

    ## get grid_mapping projection key
    try:
        #pkey = getattr(nc,'projection_key')
        pkey = nc.getncattr('projection_key')
    except:
        pkey = None

    ## write data
    if dataset in nc.variables.keys():
        ## dataset already in NC file
        ## update existing dataset attributes
        if dataset_attributes is not None:
            for att in dataset_attributes.keys():
                if att in ['_FillValue']: continue
                nc.variables[dataset].setncattr(att, dataset_attributes[att])
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
        if dataset in ['lat', 'lon']:
            netcdf_least_significant_digit = None
        else:
            netcdf_least_significant_digit = None if netcdf_compression_least_significant_digit is None else 1 * netcdf_compression_least_significant_digit

        ## select dataset type, default type or target type if discretising
        dtype = data.dtype if (not discretise) else pdisc['target_type']

        ## create variable
        var = nc.createVariable(dataset,dtype,('y','x'), fill_value = fillvalue, chunksizes = chunksizes,
                                zlib = netcdf_compression, complevel = netcdf_compression_level,
                                least_significant_digit = netcdf_least_significant_digit)

        ## add discretisation
        if discretise:
            var.scale_factor = pdisc['scale_factor']
            var.add_offset = pdisc['add_offset']
            var.set_auto_scale(True)

        if wavelength is not None: var.setncattr('wavelength', float(wavelength))
        ## set attributes
        if dataset_attributes is not None:
            for att in dataset_attributes.keys():
                if att in ['_FillValue', 'scale_factor', 'add_offset']: continue
                var.setncattr(att, dataset_attributes[att])

        ## add grid mapping key if there is projection set
        if pkey is not None: var.setncattr('grid_mapping', pkey)

        if offset is None:
            if data.dtype in (np.float32, np.float64): var[:] = np.nan
            var[:] = data
        else:
            if data.dtype in (np.float32, np.float64): var[:] = np.nan
            var[offset[1]:offset[1]+dims[0],offset[0]:offset[0]+dims[1]] = data
    if keep is not True: data = None

    ## return NC dataset
    if return_nc: return(nc)

    ## close netcdf file
    if open_file:
        nc.close()
        nc=None
