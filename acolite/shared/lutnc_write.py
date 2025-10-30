## lutnc_write
## generic netcdf lut writing
##
## written by Quinten Vanhellemont, RBINS
## 2022-03-23
## modifications: 2025-10-30 (QV) added directory check

def lutnc_write(lutnc, lut, meta_in, dims = None, override = True, compression = True, complevel = 1):
    import os
    from netCDF4 import Dataset
    import numpy as np
    import acolite as ac

    if os.path.exists(lutnc) & override:
        os.remove(lutnc)

    ## create output directory
    if not os.path.exists(os.path.dirname(lutnc)):
        os.makedirs(os.path.dirname(lutnc))

    ## update dims
    meta = {k: meta_in[k] for k in meta_in}
    if ('dims' not in meta) & (dims != None): meta['dims'] = dims
    if 'dims' in meta: dims = meta['dims']

    ## open NetCDF file in write mode
    nc = Dataset(lutnc, 'w', format='NETCDF4_CLASSIC')

    ## set attributes
    for i in meta:
        attdata=meta[i]
        if isinstance(attdata,list):
            if isinstance(attdata[0],str):
                attdata=','.join(attdata)
        setattr(nc, i, attdata)

    ## create dimensions
    for d in dims:
        nc.createDimension(d, len(meta[d]))

    ## write lut
    if type(lut) is dict:
        for k in lut.keys():
            var = nc.createVariable(k,np.float32,(dims), zlib=compression, complevel=complevel)
            nc.variables[k][:] = lut[k].astype(np.float32)
    else:
        var = nc.createVariable('lut',np.float32,(dims), zlib=compression, complevel=complevel)
        var[:] = lut.astype(np.float32)
    nc.close()


    return(lutnc)
