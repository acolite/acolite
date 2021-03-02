## import sky reflectance lut
## QV 2020-03-17
## last updates: 2020-07-15 (QV) added sensor resampling
##               2020-07-25 (QV) new default LUT
##               2021-01-14 (QV) changed input of RSKY LUT
##               2021-01-19 (QV) added support for bzipped LUTs
##               2021-01-26 (QV) new default LUT, flip the raa from OSOAA
##               2021-02-24 (QV) renamed from rsky_read_lut
##               2021-03-01 (QV) removed separate luts for wind speed

def import_rsky_lut(model, lutbase='ACOLITE-RSKY-202102-82W', sensor=None, override=False):
    import os
    import numpy as np
    import scipy.interpolate
    from netCDF4 import Dataset
    import acolite as ac

    if True:
            lutdir = '{}/LUT/{}/'.format(ac.config['data_dir'], '-'.join(lutbase.split('-')[1:3]))
            lutnc = '{}/{}-MOD{}.nc'.format(lutdir, lutbase, model)

            ## base lut
            if sensor is None:
                ## extract bz2 files
                unzipped = False
                lutncbz2 = '{}.bz2'.format(lutnc)
                if (not os.path.isfile(lutnc)) & (os.path.isfile(lutncbz2)):
                    import bz2, shutil
                    with bz2.BZ2File(lutncbz2) as fi, open(lutnc,"wb") as fo:
                        shutil.copyfileobj(fi,fo)
                    unzipped = True
                ## end extract bz2 files

                nc = Dataset(lutnc)
                meta = {}
                for attr in nc.ncattrs():
                    attdata = getattr(nc,attr)
                    if isinstance(attdata,str): attdata = attdata.split(',')
                    meta[attr]=attdata
                lut = nc.variables['lut'][:]
                lut = np.flip(lut, axis=1) ## flip raa
                nc.close()
                if unzipped: os.remove(lutnc) ## clear unzipped LUT

                dim = [meta['wave'], meta['azi'], meta['thv'], meta['ths'], meta['tau']]
                #if 'press' in meta:
                #    if type(meta['press']) is np.ndarray:
                #        dim = [meta['press'], meta['wave'], meta['azi'], meta['thv'], meta['ths'], meta['tau']]
                if 'wind' in meta:
                    if type(meta['wind']) is np.ndarray:
                        dim = [meta['wave'], meta['azi'], meta['thv'], meta['ths'], meta['wind'], meta['tau']]

                rgi = scipy.interpolate.RegularGridInterpolator(dim, lut, bounds_error=False, fill_value=np.nan)
                return(lut, meta, dim, rgi)
            ## sensor specific lut
            else:
                lutnc_s = '{}/{}/{}-MOD{}_{}.nc'.format(lutdir, sensor, lutbase, model, sensor)

                ## get sensor RSR
                pp_path = ac.config['data_dir']
                lutdir=pp_path+'/LUT/'
                rsr_file = pp_path+'/RSR/'+sensor+'.txt'
                rsr, rsr_bands = ac.shared.rsr_read(file=rsr_file)

                ## make new sensor lutfile
                if (override) & (os.path.isfile(lutnc_s)): os.remove(lutnc_s)

                if not os.path.isfile(lutnc_s):
                    print('Resampling RSKY LUT {} to {}'.format(os.path.basename(lutnc), sensor))
                    if not os.path.exists(os.path.dirname(lutnc_s)): os.makedirs(os.path.dirname(lutnc_s))
                    ## read lut
                    lut, meta, dim, rgi = ac.aerlut.import_rsky_lut(model, lutbase=lutbase)
                    ## resample to bands
                    lut_sensor = {}
                    for band in rsr_bands:
                        lut_sensor[band] = ac.shared.rsr_convolute_nd(lut, meta['wave'],rsr[band]['response'], rsr[band]['wave'], axis=0)
                    #return(lut_sensor, meta, dim)
                    ## save to new file
                    from netCDF4 import Dataset
                    nc = Dataset(lutnc_s, 'w', format='NETCDF4_CLASSIC')
                    for i in meta.keys():
                        attdata=meta[i]
                        if isinstance(attdata,list):
                            if isinstance(attdata[0],str):
                                attdata=','.join(attdata)
                        setattr(nc, i, attdata)
                    nc.createDimension('azi', len(dim[1]))
                    nc.createDimension('thv', len(dim[2]))
                    nc.createDimension('ths', len(dim[3]))
                    nc.createDimension('wind', len(dim[4]))
                    nc.createDimension('tau', len(dim[5]))
                    ## write LUT
                    for band in rsr_bands:
                        var = nc.createVariable(band,np.float32,('azi','thv','ths','wind', 'tau'))
                        nc.variables[band][:] = lut_sensor[band].astype(np.float32)
                    nc.close()
                ## end resample lut

                ## lutfile was already resampled
                if os.path.isfile(lutnc_s):
                    nc = Dataset(lutnc_s)
                    meta = {}
                    for attr in nc.ncattrs():
                        attdata = getattr(nc,attr)
                        if isinstance(attdata,str): attdata = attdata.split(',')
                        meta[attr]=attdata
                    lut_sensor = {}
                    for band in rsr_bands:
                        lut_sensor[band] = nc.variables[band][:]
                    nc.close()
                    dim = [meta['azi'], meta['thv'], meta['ths'], meta['wind'], meta['tau']]
                    rgi_sensor = {}
                    for band in rsr_bands:
                        rgi_sensor[band] = scipy.interpolate.RegularGridInterpolator(dim, lut_sensor[band], bounds_error=False, fill_value=np.nan)
                    return(lut_sensor, meta, dim, rgi_sensor)
