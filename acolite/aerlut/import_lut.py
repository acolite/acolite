## def import_lut
## imports LUT made with 6SV and converts to NetCDF
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07-05
## modifications:   2020-07-14 (QV)
##                  2021-01-16 (QV) added support for bz2 compressed luts
##                  2021-02-24 (QV) removed obsolete code
##                  2021-02-25 (QV) changed position of lut files (removed lutid directory), added removal of unzipped file
##                  2021-03-02 (QV) integrated sensor specific LUTs
##                  2021-05-31 (QV) added remote lut retrieval
##                  2021-06-08 (QV) added lut par subsetting
##                  2021-07-20 (QV) added retrieval of generic LUTs
##                  2021-10-22 (QV) compute ttot if not in LUT

def import_lut(lutid, lutdir, lut_par = ['utott', 'dtott', 'astot', 'ttot', 'romix'],
               override = False, sensor = None, get_remote = True,
               remote_base = 'https://raw.githubusercontent.com/acolite/acolite_luts/main'):

    import os, sys
    import numpy as np
    import acolite as ac

    lutnc=lutdir+'/'+lutid+'.nc'
    lut = None

    ## generic LUT
    if sensor is None:
        ## extract bz2 files
        unzipped = False
        lutncbz2 = '{}.bz2'.format(lutnc)

        ## try downloading LUT from GitHub
        if (not os.path.isfile(lutnc)) & (not os.path.isfile(lutncbz2)) & (get_remote):
            remote_lut = '{}/{}/{}'.format(remote_base, '-'.join(lutid.split('-')[0:3]), os.path.basename(lutncbz2))
            try:
                print('Getting remote LUT {}'.format(remote_lut))
                ac.shared.download_file(remote_lut, lutncbz2)
            except:
                print('Could not download remote lut {} to {}'.format(remote_lut, lutncbz2))
                if os.path.exists(lutncbz2): os.remove(lutncbz2)

        ## extract bz LUT
        if (not os.path.isfile(lutnc)) & (os.path.isfile(lutncbz2)):
            import bz2, shutil
            with bz2.BZ2File(lutncbz2) as fi, open(lutnc,"wb") as fo:
                shutil.copyfileobj(fi,fo)
            unzipped = True
        ## end extract bz2 files

        ## read dataset from NetCDF
        try:
            lut, meta = ac.shared.lutnc_import(lutnc)
        except:
            print(sys.exc_info()[0])
            print('Failed to open LUT data from NetCDF (id='+lutid+')')

        if unzipped: os.remove(lutnc) ## clear unzipped LUT

        if lut is None:
            print('Could not import LUT {} from {}'.format(lutid, lutdir))
            return()

        ## subset LUTs
        if lut_par is not None:
            lut_sub_idx = []
            lut_sub_par = []
            for i, ik in enumerate(lut_par):
                for j, jk in enumerate(meta['par']):
                    if (ik == jk):
                        lut_sub_idx.append(j)
                        lut_sub_par.append(jk)
            ## add ttot if not in NetCDF
            if ('ttot' in lut_par) & ('ttot' not in meta['par']):
                for j, jk in enumerate(meta['par']):
                    if 'tray' == jk: tri = j
                    if 'taer' == jk: tai = j
                lut = np.append(lut, lut[[tri], :,:,:,:,:,:]+lut[[tai], :,:,:,:,:,:], axis=0)
                lut_sub_par.append('ttot')
                lut_sub_idx.append(lut.shape[0]-1)

            ## add romix+rsurf if not in NetCDF
            if ('romix+rsurf' in lut_par) & ('romix+rsurf' not in meta['par']) & \
               (('romix' in meta['par']) & ('rsurf' in meta['par'])):
                for j, jk in enumerate(meta['par']):
                    if 'romix' == jk: pi = j
                    if 'rsurf' == jk: pj = j
                lut = np.append(lut, lut[[pi], :,:,:,:,:,:]+lut[[pj], :,:,:,:,:,:], axis=0)
                lut_sub_par.append('romix+rsurf')
                lut_sub_idx.append(lut.shape[0]-1)

            ## subset to requested parameters
            meta['par'] = lut_sub_par
            lut = lut[lut_sub_idx,:,:,:,:,:,:]

        ## for the  and Continental and Urban models (1,3)
        ## romix nans were retrieved for wavelengths > 2 micron and aot == 0.001
        ## 500mb for C+U and 1013/1100 for U
        ## if any nans set then to 0
        sub = np.where(np.isnan(lut))
        lut[sub] = 0
        return(lut, meta)

    ## sensor specific LUT
    else:
        ## sensor LUT NetCDF is stored here
        lutnc_s='{}/{}/{}_{}.nc'.format(lutdir,sensor,lutid,sensor)
        if not os.path.exists(os.path.dirname(lutnc_s)): os.makedirs(os.path.dirname(lutnc_s))
        if (os.path.isfile(lutnc_s)) & (override): os.remove(lutnc_s)

        if (not os.path.isfile(lutnc_s)) | (override):
            ## try downloading LUT from GitHub
            if (get_remote):
                slut = '{}_{}'.format(lutid, sensor)
                remote_lut = '{}/{}/{}/{}.nc'.format(remote_base, '-'.join(lutid.split('-')[0:3]), sensor, slut)
                try:
                    print('Getting remote LUT {}'.format(remote_lut))
                    ac.shared.download_file(remote_lut, lutnc_s)
                    print('Testing LUT {}'.format(lutnc_s))
                    lut, meta = ac.shared.lutnc_import(lutnc_s) # test LUT
                except:
                    print('Could not download remote lut {} to {}'.format(remote_lut, lutnc_s))
                    if os.path.exists(lutnc_s): os.remove(lutnc_s)

            ## otherwise to local resampling
            if (not os.path.isfile(lutnc_s)):
                print('Resampling LUT {} to sensor {}'.format(lutid, sensor))
                rsrd = ac.shared.rsr_dict(sensor=sensor)
                rsr, rsr_bands = rsrd[sensor]['rsr'], rsrd[sensor]['rsr_bands']

                ## read LUT
                lut, meta = ac.aerlut.import_lut(lutid,lutdir, lut_par=None) ## add None so all pars are loaded when resampling
                lut_dims = lut.shape

                ## new ndim convolution
                lut_sensor = {}
                for band in rsr_bands:
                    lut_sensor[band] = ac.shared.rsr_convolute_nd(lut, meta['wave'],rsr[band]['response'], rsr[band]['wave'], axis=1)

                ## write nc file
                try:
                    if os.path.isfile(lutnc_s) is False:
                        from netCDF4 import Dataset
                        nc = Dataset(lutnc_s, 'w', format='NETCDF4_CLASSIC')
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
                            var = nc.createVariable(band,np.float32,('par','azi','thv','ths','wnd','tau'))
                            nc.variables[band][:] = lut_sensor[band].astype(np.float32)
                        nc.close()
                        nc = None
                        arr = None
                        meta = None
                except:
                    if os.path.isfile(lutnc_s): os.remove(lutnc_s)
                    print(sys.exc_info()[0])
                    print('Failed to write LUT data to NetCDF (id='+lutid+')')

        ## read dataset from NetCDF
        if os.path.isfile(lutnc_s):
            try:
                lut_sensor, meta = ac.shared.lutnc_import(lutnc_s)
            except:
                print(sys.exc_info()[0])
                print('Failed to open LUT data from NetCDF (id='+lutid+')')

        ## subset LUTs
        if lut_par is not None:
            lut_sub_idx = []
            lut_sub_par = []
            for i, ik in enumerate(lut_par):
                for j, jk in enumerate(meta['par']):
                    if (ik == jk):
                        lut_sub_idx.append(j)
                        lut_sub_par.append(jk)
            ## add ttot if not in NetCDF
            if ('ttot' in lut_par) & ('ttot' not in meta['par']):
                for j, jk in enumerate(meta['par']):
                    if 'tray' == jk: tri = j
                    if 'taer' == jk: tai = j
                for dataset in lut_sensor:
                    lut_sensor[dataset] = np.append(lut_sensor[dataset], \
                        lut_sensor[dataset][[tri],:,:,:,:,:]+lut_sensor[dataset][[tai],:,:,:,:,:], axis=0)
                lut_sub_par.append('ttot')
                lut_sub_idx.append(lut_sensor[dataset].shape[0]-1)

            ## add romix+rsurf if not in NetCDF
            if ('romix+rsurf' in lut_par) & ('romix+rsurf' not in meta['par']) &\
               (('romix' in meta['par']) & ('rsurf' in meta['par'])):
                for j, jk in enumerate(meta['par']):
                    if 'romix' == jk: pi = j
                    if 'rsurf' == jk: pj = j
                for dataset in lut_sensor:
                    lut_sensor[dataset] = np.append(lut_sensor[dataset], \
                        lut_sensor[dataset][[pi],:,:,:,:,:]+lut_sensor[dataset][[pj],:,:,:,:,:], axis=0)
                lut_sub_par.append('romix+rsurf')
                lut_sub_idx.append(lut_sensor[dataset].shape[0]-1)

            meta['par'] = lut_sub_par
            for dataset in lut_sensor:
                lut_sensor[dataset] = lut_sensor[dataset][lut_sub_idx,:,:,:,:,:]

        return(lut_sensor, meta)
