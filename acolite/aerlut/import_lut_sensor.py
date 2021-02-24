## def import_sensor_lut
## imports LUT, interpolates to hyperspectral and convolutes to sensor bands
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07-13
## modifications: 2017-01-23 (QV) added generic lutdir/rsrfile locations
##                2018-06-07 (QV) reordered lutdir and rsr_file check
##                2018-07-18 (QV) changed acolite import name
##                2020-07-14 (QV) cleaned up, new interpolation, and SWIR nans are set to 0
##                2021-02-24 (QV) cleaned up, renamed from get_sensor_lut

def import_lut_sensor(sensor, rsr_file, lutid, override=0, lutdir=None):
    import os, sys
    import numpy as np
    import acolite as ac

    ## get sensor RSR
    if lutdir is None: lutdir=ac.config['data_dir']+'/LUT/'

    ## sensor LUT NetCDF is stored here
    lutnc=lutdir+'/'+lutid+'/'+lutid+'_'+sensor+'.nc'
    if (os.path.isfile(lutnc)) & (override == 1): os.remove(lutnc)

    ## generate sensor LUT NetCDF if needed
    if (os.path.isfile(lutnc) is False):
        print('Resampling LUT {} to sensor {}'.format(lutid, sensor))
        if rsr_file is None: rsr_file = ac.config['data_dir']+'/RSR/'+sensor+'.txt'
        rsr, rsr_bands = ac.shared.rsr_read(file=rsr_file)

        ## read LUT
        lut, meta = ac.aerlut.import_lut(lutid,lutdir, override=0)
        lut_dims = lut.shape

        ## new ndim convolution
        lut_sensor = {}
        for band in rsr_bands:
            lut_sensor[band] = ac.shared.rsr_convolute_nd(lut, meta['wave'],rsr[band]['response'], rsr[band]['wave'], axis=1)

        ## write nc file
        try:
            if os.path.isfile(lutnc) is False:
                from netCDF4 import Dataset
                nc = Dataset(lutnc, 'w', format='NETCDF4_CLASSIC')
                ## write metadata
                for i in meta:
                    attdata=meta[i]
                    if isinstance(attdata,list):
                        if isinstance(attdata[0],str):
                            attdata=','.join(attdata)
                    setattr(nc, i, attdata)
                ## set up LUT dimension
                nc.createDimension('par', lut_dims[0])
                #nc.createDimension('wave', lut_dims[1]) # not used here
                nc.createDimension('azi', lut_dims[2])
                nc.createDimension('thv', lut_dims[3])
                nc.createDimension('ths', lut_dims[4])
                nc.createDimension('wnd', lut_dims[5])
                nc.createDimension('tau', lut_dims[6])
                ## write LUT
                for band in lut_sensor.keys():
                    var = nc.createVariable(band,float,('par','azi','thv','ths','wnd','tau'))
                    nc.variables[band][:] = lut_sensor[band]
                nc.close()
                nc = None
                arr = None
                meta = None
        except:
            if os.path.isfile(lutnc): os.remove(lutnc)
            print(sys.exc_info()[0])
            print('Failed to write LUT data to NetCDF (id='+lutid+')')

    ## read dataset from NetCDF
    if os.path.isfile(lutnc):
        try:
            if os.path.isfile(lutnc):
                from netCDF4 import Dataset
                nc = Dataset(lutnc)
                ## read in metadata
                meta=dict()
                for attr in nc.ncattrs():
                    attdata = getattr(nc,attr)
                    if isinstance(attdata,str): attdata = attdata.split(',')
                    meta[attr]=attdata
                ## read in LUT
                lut_sensor = dict()
                datasets = list(nc.variables.keys())
                for dataset in datasets:
                    lut_sensor[dataset] = nc.variables[dataset][:]
                nc.close()
                nc = None
        except:
            print(sys.exc_info()[0])
            print('Failed to open LUT data from NetCDF (id='+lutid+')')

    return(lut_sensor, meta)
