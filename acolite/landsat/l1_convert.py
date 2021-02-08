## def l1_convert
## converts landsat bundle data to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-02-05
## modifications: 2021-02-07 (QV) added support for tile merging and extending the size of the output file to the requested limit

def l1_convert(inputfile, output=None,
                limit=None, sub=None,
                output_pan = True,
                output_pan_ms = True,
                output_thermal = True,

                output_geometry = True,
                output_geolocation = True,
                output_xy = False,
                usgs_reflectance = True,

                percentiles_compute = True,
                percentiles = (0,1,5,10,25,50,75,90,95,99,100),

                merge_tiles = False,
                extend_region = False,

                verbosity = 0, vname = ''):

    import os, glob, dateutil, time
    import acolite as ac
    import scipy.ndimage
    import numpy as np
    t0 = time.time()

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ## check if merging settings make sense
    if (limit is None) & (merge_tiles):
        if verbosity > 0: print("Merging tiles not supported without ROI limit")
        merge_tiles = False
    if merge_tiles: extend_region = True

    new = True
    new_pan = True
    ofile = None
    ofiles = []

    for bundle in inputfile:
        if output is None: output = os.path.dirname(bundle)
        if verbosity > 1: print('Starting conversion of {}'.format(bundle))

        mtl = glob.glob('{}/{}'.format(bundle, '*MTL.txt'))
        if len(mtl) == 0:
            if verbosity > 0: print('No metadata file found for {}'.format(bundle))
            return(1)
        else:
            mtl = mtl[0]

        ## read landsat metadata and check files
        if verbosity > 1: print('Importing metadata from {}'.format(bundle))
        meta = ac.landsat.metadata_read(mtl)
        fmeta = ac.landsat.metadata_bands(bundle, meta)

        ## get relevant data from meta
        if 'PRODUCT_CONTENTS' in meta: ## COLL2
            pk = 'IMAGE_ATTRIBUTES'
            ik = 'IMAGE_ATTRIBUTES'
            rk = 'PROJECTION_ATTRIBUTES'
        elif 'PRODUCT_METADATA' in meta: ## COLL1
            pk = 'PRODUCT_METADATA'
            ik = 'IMAGE_ATTRIBUTES'
            rk = 'PRODUCT_METADATA'
        spacecraft_id = meta[pk]['SPACECRAFT_ID']
        sensor_id = meta[pk]['SENSOR_ID']
        path = meta[pk]['WRS_PATH']
        row = meta[pk]['WRS_ROW']
        isodate = meta[pk]['DATE_ACQUIRED']+'T'+meta[pk]['SCENE_CENTER_TIME']
        global_dims = int(meta[rk]['REFLECTIVE_LINES']), int(meta[rk]['REFLECTIVE_SAMPLES'])

        ## some hard coded info
        sat = 'L{}'.format(spacecraft_id[-1])
        if sensor_id in ['TM']:# Landsat 5
            sen = 'TM'
            pan_bands = []
            thermal_bands = ['6']
        elif sensor_id in ['ETM']:# Landsat 7
            sen = 'ETM'
            pan_bands = ['8']
            thermal_bands = ['6_VCID_1', '6_VCID_2']
        elif sensor_id in ['OLI', 'OLI_TIRS']:# Landsat 8
            sen = 'OLI'
            pan_bands = ['8']
            thermal_bands = ['10', '11']
        else:
            print(spacecraft_id, sensor_id)
            print('Not configured')
            return(1)

        sensor = '{}_{}'.format(sat,sen)

        ## scene average geometry
        vza = 0
        sza = 90-float(meta[ik]['SUN_ELEVATION'])
        raa = float(meta[ik]['SUN_AZIMUTH'])
        se_distance = float(meta[ik]['EARTH_SUN_DISTANCE'])

        ## read rsr
        rsrf = ac.path+'/data/RSR/{}.txt'.format(sensor)
        rsr, rsr_bands = ac.shared.rsr_read(rsrf)
        waves = np.arange(250, 2500)/1000
        waves_mu = ac.shared.rsr_convolute_dict(waves, waves, rsr)
        waves_names = {'{}'.format(b):'{:.0f}'.format(waves_mu[b]*1000) for b in waves_mu}

        ## get F0 - not stricty necessary if using USGS reflectance
        f0 = ac.shared.f0_get()
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data'])*10, rsr)

        ## make global attributes for L1R NetCDF
        gatts = {'sensor':sensor, 'isodate':isodate, 'global_dims':global_dims,
                 'sza':sza, 'vza':vza, 'raa':raa, 'se_distance': se_distance,
                 'mus': np.cos(sza*(np.pi/180.)), 'wrs_path': path, 'wrs_row': row}
        if merge_tiles:
            gatts['tile_code'] = 'merged'
        else:
            gatts['tile_code'] = '{}{}'.format(gatts['wrs_path'].zfill(3),gatts['wrs_row'].zfill(3))

        stime = dateutil.parser.parse(gatts['isodate'])
        oname = '{}_{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'), gatts['tile_code'])
        if vname != '': oname+='_{}'.format(vname)

        ## output file information
        if (merge_tiles is False) | (ofile is None):
            ofile = '{}/{}_L1R.nc'.format(output, oname)
            gatts['oname'] = oname
            gatts['ofile'] = ofile
        elif (merge_tiles) & (ofile is None):
            ofile = '{}/{}_L1R.nc'.format(output, oname)
            gatts['oname'] = oname
            gatts['ofile'] = ofile

        print(ofile)

        ## add band info to gatts
        for b in rsr_bands:
            gatts['{}_wave'.format(b)] = waves_mu[b]*1000
            gatts['{}_name'.format(b)] = waves_names[b]
            gatts['{}_f0'.format(b)] = f0_b[b]
            if b in fmeta:
                fmeta[b]['f0'] = f0_b[b]
                fmeta[b]['se_distance'] = gatts['se_distance']

        ## get scene projection and extent
        dct = ac.landsat.projection(meta)
        gatts['xrange'] = dct['xrange']
        gatts['yrange'] = dct['yrange']
        gatts['proj4_string'] = dct['proj4_string']
        gatts['pixel_size'] = dct['pixel_size']
        if 'zone' in dct: gatts['zone'] = dct['zone']

        if (sub is None) & (limit is not None):
            dct_sub = ac.shared.projection_sub(dct, limit, four_corners=True)
            if dct_sub['out_lon']:
                if verbosity > 1: print('Longitude limits outside {}'.format(bundle))
                return(1)
            if dct_sub['out_lat']:
                if verbosity > 1: print('Latitude limits outside {}'.format(bundle))
                return(1)
            sub = dct_sub['sub']
        else:
            if extend_region:
                print("Can't extend region if no ROI limits given")
                extend_region = False

        ## get pan dimensions and subset
        nc_dim = None
        nc_off = None
        #pan_padding = None
        if sub is None:
            pan_ms_dims = global_dims
            pan_dims = global_dims[0]*2, global_dims[1]*2
            #if 'PANCHROMATIC_LINES' in meta[rk]:
            #    pan_dims_ = int(meta[rk]['PANCHROMATIC_LINES']), int(meta[rk]['PANCHROMATIC_SAMPLES'])
            #    pan_padding = pan_dims[0]-pan_dims_[0], pan_dims[1]-pan_dims_[1]
            sub_pan = None
            gatts['pan_dims'] = pan_dims
            gatts['global_dims'] = dct['dimensions']
        else:
            pan_ms_dims = sub[3], sub[2]
            pan_dims = sub[3]*2, sub[2]*2
            sub_pan = [s*2 for s in sub]
            gatts['sub'] = sub
            gatts['global_dims'] = dct_sub['ydim'], dct_sub['xdim']
            gatts['scene_dims'] = dct['ydim'], dct['xdim']
            gatts['limit'] = limit
            gatts['pan_dims'] = pan_dims
            gatts['pan_sub'] = sub_pan

            if extend_region:
                nc_dim = dct_sub['region']['ydim'], dct_sub['region']['xdim']
                nc_off_yx = (int((dct_sub['yrange'][0]-dct_sub['region']['yrange'][0])/dct_sub['pixel_size'][1]),
                             int((dct_sub['xrange'][0]-dct_sub['region']['xrange'][0])/dct_sub['pixel_size'][0]))
                nc_off = (int((dct_sub['xrange'][0]-dct_sub['region']['xrange'][0])/dct_sub['pixel_size'][0]),
                          int((dct_sub['yrange'][0]-dct_sub['region']['yrange'][0])/dct_sub['pixel_size'][1]))
                gatts['global_dims'] = nc_dim

        ## new file if not merging
        if (merge_tiles is False):
            new = True
            new_pan = True

        ## start the conversion
        ## write geometry
        if ('VAA' in fmeta) & ('SAA' in fmeta) & ('VZA' in fmeta) & ('SZA' in fmeta):
            if verbosity > 1: print('Reading per pixel geometry')
            sza = ac.landsat.read_band(fmeta['SZA']['FILE'], sub=sub).astype(np.float32)/100
            mus = np.cos(sza*(np.pi/180.)) ## per pixel cos sun zenith
            if (output_geometry):
                saa = ac.landsat.read_band(fmeta['SAA']['FILE'], sub=sub).astype(np.float32)/100
                vza = ac.landsat.read_band(fmeta['VZA']['FILE'], sub=sub).astype(np.float32)/100
                vaa = ac.landsat.read_band(fmeta['VAA']['FILE'], sub=sub).astype(np.float32)/100
                mask = (vaa == 0) * (vza == 0) * (saa == 0) * (sza == 0)
                vza[mask] = np.nan
                sza[mask] = np.nan
                #vaa[mask] = np.nan
                #saa[mask] = np.nan

                raa = (saa-vaa)
                raa[raa>180]-=180
                raa[mask] = np.nan
                vaa = None
                saa = None
                mask = None
                ac.output.nc_write(ofile, 'raa', raa, global_dims=nc_dim, offset=nc_off, replace_nan=True, attributes=gatts, new=new)
                if verbosity > 1: print('Wrote raa')
                new = False
                ac.output.nc_write(ofile, 'vza', vza, global_dims=nc_dim, offset=nc_off, replace_nan=True)
                if verbosity > 1: print('Wrote vza')
                ac.output.nc_write(ofile, 'sza', sza, global_dims=nc_dim, offset=nc_off, replace_nan=True)
                if verbosity > 1: print('Wrote sza')

                #ac.output.nc_write(ofile, 'saa', saa, global_dims=nc_dim, offset=nc_off, replace_nan=True)
                #if verbosity > 1: print('Wrote saa')
                #ac.output.nc_write(ofile, 'vaa', vaa, global_dims=nc_dim, offset=nc_off, replace_nan=True)
                #if verbosity > 1: print('Wrote vaa')
                sza = None
                vza = None
        else:
            mus = np.asarray(gatts['mus'])  ## average cos sun zenith
            #mus.shape+=(1,1)

        ## write lat/lon
        if (output_geolocation):
            if os.path.exists(ofile):
                datasets = ac.shared.nc_datasets(ofile)
            else:
                datasets = []
            if ('lat' not in datasets) or ('lon' not in datasets):
                if verbosity > 1: print('Writing geolocation lon/lat')
                if extend_region:
                    lon, lat = ac.shared.projection_geo(dct_sub['region'])
                else:
                    lon, lat = ac.shared.projection_geo(dct if sub is None else dct_sub)
                ac.output.nc_write(ofile, 'lon', lon, global_dims=nc_dim, attributes=gatts, new=new, double=True)
                if verbosity > 1: print('Wrote lon')
                ac.output.nc_write(ofile, 'lat', lat, global_dims=nc_dim, double=True)
                if verbosity > 1: print('Wrote lat')
                new=False

        ## write x/y
        if (output_xy):
            if os.path.exists(ofile):
                datasets = ac.shared.nc_datasets(ofile)
            else:
                datasets = []
            if ('x' not in datasets) or ('y' not in datasets):
                if verbosity > 1: print('Writing geolocation x/y')
                if extend_region:
                    x, y = ac.shared.projection_geo(dct_sub['region'], xy=True)
                else:
                    x, y = ac.shared.projection_geo(dct if sub is None else dct_sub, xy=True)
                #x, y = ac.shared.projection_geo(dct if sub is None else dct_sub, xy=True)
                ac.output.nc_write(ofile, 'x', x, global_dims=nc_dim, attributes=gatts, new=new)
                if verbosity > 1: print('Wrote x')
                ac.output.nc_write(ofile, 'y', y, global_dims=nc_dim)
                if verbosity > 1: print('Wrote y')
                new=False

        ## write TOA bands
        if verbosity > 1: print('Converting bands')
        for b in fmeta:
            if '.TIF' not in fmeta[b]['FILE']: continue
            if b in ['PIXEL', 'RADSAT']: continue
            if os.path.exists(fmeta[b]['FILE']):
                if b in waves_names:
                    pan = False
                    if b in pan_bands:
                        if (not output_pan) & (not output_pan_ms): continue
                        pan = True
                        if len(np.atleast_1d(mus))>1:
                            mus_pan = scipy.ndimage.zoom(mus, zoom=2, order=1)
                        # this is now done in read_toa? to also support crops that go over the border
                        #    ## subset if we are "missing" a pan pixel
                        #    if pan_padding is not None:
                        #        mus_pan = mus_pan[0:pan_dims[0]-pan_padding[0],
                        #                          0:pan_dims[1]-pan_padding[1]]

                        else:
                            mus_pan = mus * 1
                        data = ac.landsat.read_toa(fmeta[b], sub=sub_pan, mus=mus_pan)
                        mus_pan = None
                    else:
                        data = ac.landsat.read_toa(fmeta[b], sub=sub, mus=mus)

                    ds = 'rhot_{}'.format(waves_names[b])
                    ds_att = {'wavelength':waves_mu[b]*1000}
                    for k in fmeta[b]: ds_att[k] = fmeta[b][k]
                    if percentiles_compute:
                        ds_att['percentiles'] = percentiles
                        ds_att['percentiles_data'] = np.nanpercentile(data, percentiles)

                    if output_pan & pan:
                        ofile_pan = ofile.replace('_L1R.nc', '_L1R_pan.nc')

                        ## add padding to pan data
                        if data.shape[0] <  pan_dims[0]:
                            data = np.vstack((data, np.zeros((pan_dims[0]-data.shape[0], data.shape[1]))))
                        elif data.shape[0] >  pan_dims[0]:
                            data = data[0:pan_dims[0], :]
                        if data.shape[1] < pan_dims[1]:
                            data = np.hstack((data, np.zeros((data.shape[0], pan_dims[1]-data.shape[1]))))
                        elif data.shape[1] > pan_dims[1]:
                            data = data[:, 0:pan_dims[1]]

                        ## write output
                        nc_dim_pan = None
                        nc_off_pan = None
                        if (nc_dim is not None): ## this means we are "extending" the output size
                            nc_dim_pan = nc_dim[0]*2, nc_dim[1]*2
                            nc_off_pan = nc_off[0]*2, nc_off[1]*2
                        ac.output.nc_write(ofile_pan, ds, data, attributes=gatts,
                                           global_dims=nc_dim_pan, offset=nc_off_pan, replace_nan=True,
                                           new=new_pan, dataset_attributes = ds_att)
                        new_pan = False
                        if verbosity > 1: print('Converting bands: Wrote {} to separate L1R_pan'.format(ds))

                        if output_pan_ms:
                            ## prepare for low res output
                            data = scipy.ndimage.zoom(data, zoom=0.5, order=1)
                        else:
                            continue

                    #if data.shape == gatts['global_dims']:
                    if True:
                        ac.output.nc_write(ofile, ds, data,
                                           global_dims=nc_dim, offset=nc_off, replace_nan=True,
                                           attributes=gatts, new=new, dataset_attributes = ds_att)
                        new = False
                        if verbosity > 1: print('Converting bands: Wrote {}'.format(ds))
                    else:
                        if verbosity > 0: print('Converting bands: Error in writing {}'.format(ds))
                else:
                    if b in thermal_bands:
                        if output_thermal:
                            ds = 'bt{}'.format(b).lower()
                            ds_att = {'band':b}
                            for k in fmeta[b]: ds_att[k] = fmeta[b][k]
                            if percentiles_compute:
                                ds_att['percentiles'] = percentiles
                                ds_att['percentiles_data'] = np.nanpercentile(data, percentiles)
                            data = ac.landsat.read_toa(fmeta[b], sub=sub)
                            ac.output.nc_write(ofile, ds, data,
                                               global_dims=nc_dim, offset=nc_off, replace_nan=True,
                                               attributes=gatts, new=new, dataset_attributes=ds_att)
                            new = False
                            if verbosity > 1: print('Converting bands: Wrote {}'.format(ds))
                    else:
                        continue

        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))

        if limit is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles)
