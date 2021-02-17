## def l1_convert
## converts sentinel .SAFE bundle data to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-02-11
## modifications:

def l1_convert(inputfile, output=None,
                limit=None, sub=None,
                s2_target_res = 10,

                output_geometry = True,
                output_geolocation = True,
                output_xy = False,

                geometry_type = 'gpt', ## 'gpt' or 'grids'
                geometry_res = 60, ## for gpt geometry
                geometry_format='GeoTIFF', ## for gpt geometry
                geometry_override = False, ## for gpt geometry

                percentiles_compute = True,
                percentiles = (0,1,5,10,25,50,75,90,95,99,100),

                merge_tiles = False,
                merge_zones = False,
                extend_region = False,

                check_sensor = True,
                check_time = True,
                max_merge_time = 600, # seconds

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
    if merge_tiles:
        merge_zones = True
        extend_region = True

    new = True
    warp_to = None

    ofile = None
    ofiles = []

    for bundle in inputfile:
        if output is None: output = os.path.dirname(bundle)
        if verbosity > 1: print('Starting conversion of {}'.format(bundle))

        safe_files = ac.sentinel2.safe_test(bundle)
        if 'granules' not in safe_files:
            print('File not recognised: {}'.format(bundle))
        if len(safe_files['granules']) > 1:
            print('Multi granule files are no longer supported.')
            continue

        granule = safe_files['granules'][0]

        if verbosity > 1: print('Importing metadata from {}'.format(granule))
        grmeta = ac.sentinel2.metadata_granule(safe_files[granule]['metadata']['path'])
        meta, band_data= ac.sentinel2.metadata_scene(safe_files['metadata']['path'])

        ## get relevant data from meta
        if meta['SPACECRAFT_NAME'] == 'Sentinel-2A':
            sensor = 'S2A_MSI'
        elif meta['SPACECRAFT_NAME'] == 'Sentinel-2B':
            sensor = 'S2B_MSI'
        else:
            print('{} not supported'.format(meta['SPACECRAFT_NAME']))
            continue
        if meta['PROCESSING_LEVEL'] != 'Level-1C':
            print('Processing level {} not supported'.format(meta['PROCESSING_LEVEL']))

        dtime = dateutil.parser.parse(grmeta['SENSING_TIME'])
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()
        dct = ac.sentinel2.projection(grmeta, s2_target_res=int(s2_target_res))
        global_dims = dct['dimensions']

        mgrs_tile = meta['PRODUCT_URI'].split('_')[-2]
        mgrs_tile = grmeta['TILE_ID'].split('_')[-2]

        ## scene average geometry
        sza = grmeta['SUN']['Mean_Zenith']
        saa = grmeta['SUN']['Mean_Azimuth']
        vza = np.nanmean(grmeta['VIEW']['Average_View_Zenith'])
        vaa = np.nanmean(grmeta['VIEW']['Average_View_Azimuth'])
        raa = np.abs(saa-vaa)
        while raa > 180: raa = abs(raa-360)

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
                 'mus': np.cos(sza*(np.pi/180.)), 'granule': granule, 'mgrs_tile': mgrs_tile}
        if merge_tiles:
            gatts['tile_code'] = 'merged'
        else:
            gatts['tile_code'] = '{}'.format(gatts['mgrs_tile'])

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

        ## check if we should merge these tiles
        if (merge_tiles) & (not new) & (os.path.exists(ofile)):
                fgatts = ac.shared.nc_gatts(ofile)
                if (check_sensor) & (fgatts['sensor'] != gatts['sensor']):
                    print('Sensors do not match, skipping {}'.format(bundle))
                    continue
                if check_time:
                    tdiff = dateutil.parser.parse(fgatts['isodate'])-dateutil.parser.parse(gatts['isodate'])
                    tdiff = abs(tdiff.days*86400 + tdiff.seconds)
                    if (tdiff > max_merge_time):
                        print('Time difference too large, skipping {}'.format(bundle))
                        continue

        ## add band info to gatts
        for b in rsr_bands:
            gatts['{}_wave'.format(b)] = waves_mu[b]*1000
            gatts['{}_name'.format(b)] = waves_names[b]
            gatts['{}_f0'.format(b)] = f0_b[b]
            #if b in fmeta:
            #    fmeta[b]['f0'] = f0_b[b]
            #    fmeta[b]['se_distance'] = gatts['se_distance']

        ## get scene projection and extent
        dct = ac.sentinel2.projection(grmeta, s2_target_res=int(s2_target_res))

        ## full scene
        gatts['scene_xrange'] = dct['xrange']
        gatts['scene_yrange'] = dct['yrange']
        gatts['scene_proj4_string'] = dct['proj4_string']
        gatts['scene_pixel_size'] = dct['pixel_size']
        gatts['scene_dims'] = dct['dimensions']
        if 'zone' in dct: gatts['scene_zone'] = dct['zone']

        ## check crop
        if (sub is None) & (limit is not None):
            dct_sub = ac.shared.projection_sub(dct, limit, four_corners=True)
            if dct_sub['out_lon']:
                if verbosity > 1: print('Longitude limits outside {}'.format(bundle))
                continue
            if dct_sub['out_lat']:
                if verbosity > 1: print('Latitude limits outside {}'.format(bundle))
                continue
            sub = dct_sub['sub']
        else:
            if extend_region:
                print("Can't extend region if no ROI limits given")
                extend_region = False

        ##
        if ((merge_tiles is False) & (merge_zones is False)): warp_to = None
        if sub is None:
            sub_pan = None
            if ((merge_zones) & (warp_to is not None)):
                if dct_prj != dct: ## target projection differs from this tile, need to set bounds
                    if dct['proj4_string'] != dct_prj['proj4_string']:
                        ## if the prj does not match, project current scene bounds to lat/lon
                        lonr, latr = dct['p'](dct['xrange'], dct['yrange'], inverse=True)
                        ## then to target projection
                        xrange_raw, yrange_raw = dct_prj['p'](lonr, (latr[1], latr[0]))
                        ## fix to nearest full pixel
                        pixel_size = dct_prj['pixel_size']
                        dct_prj['xrange'] = [xrange_raw[0] - (xrange_raw[0] % pixel_size[0]), xrange_raw[1]+pixel_size[0]-(xrange_raw[1] % pixel_size[0])]
                        dct_prj['yrange'] = [yrange_raw[1]+pixel_size[1]-(yrange_raw[1] % pixel_size[1]), yrange_raw[0] - (yrange_raw[0] % pixel_size[1])]
                        ## need to add new dimensions
                        dct_prj['xdim'] = int((dct_prj['xrange'][1]-dct_prj['xrange'][0])/pixel_size[0])+1
                        dct_prj['ydim'] = int((dct_prj['yrange'][1]-dct_prj['yrange'][0])/pixel_size[1])+1
                        dct_prj['dimensions'] = [dct_prj['xdim'], dct_prj['ydim']]
                    else:
                        ## if the projection matches just use the current scene projection
                        dct_prj = {k:dct[k] for k in dct}
            elif (warp_to is None):
                dct_prj = {k:dct[k] for k in dct}
        else:
            gatts['sub'] = sub
            gatts['limit'] = limit

            ## get the target NetCDF dimensions and dataset offset
            if (warp_to is None):
                if (extend_region): ## include part of the roi not covered by the scene
                    dct_prj = {k:dct_sub['region'][k] for k in dct_sub['region']}
                else: ## just include roi that is covered by the scene
                    dct_prj = {k:dct_sub[k] for k in dct_sub}
        ## end cropped

        pkeys = ['xrange', 'yrange', 'proj4_string', 'pixel_size', 'zone']
        for k in pkeys:
            if k in dct_prj: gatts[k] = dct_prj[k]

        ## else use the projection info in dct_prj
        ## from landsat processing
        #xyr = [min(dct_prj['xrange'])-dct_prj['pixel_size'][0]/2,min(dct_prj['yrange']),
        #       max(dct_prj['xrange']),max(dct_prj['yrange'])-dct_prj['pixel_size'][1]/2,
        #       dct_prj['proj4_string']]

        ## eval - this worked with another half pixel offset in computing lat lon
        xyr = [min(dct_prj['xrange']),
               min(dct_prj['yrange'])+dct_prj['pixel_size'][1],
               max(dct_prj['xrange'])+dct_prj['pixel_size'][0],
               max(dct_prj['yrange']),
               dct_prj['proj4_string']]
        ## eval -
        if False:
            xyr = [min(dct_prj['xrange']),
                   min(dct_prj['yrange'])+dct_prj['pixel_size'][1]/2,
                   max(dct_prj['xrange'])+dct_prj['pixel_size'][0]/2,
                   max(dct_prj['yrange']),
                   dct_prj['proj4_string']]
            xyr = [min(dct_prj['xrange']),
                   min(dct_prj['yrange']),
                   max(dct_prj['xrange']),
                   max(dct_prj['yrange']),
                   dct_prj['proj4_string']]

        ## warp settings for read_band
        res_method = 'average'
        warp_to = (dct_prj['proj4_string'], xyr, dct_prj['pixel_size'][0],dct_prj['pixel_size'][1], res_method)

        ## store scene and output dimensions
        gatts['scene_dims'] = dct['ydim'], dct['xdim']
        gatts['global_dims'] = dct_prj['dimensions']

        ## new file for every bundle if not merging
        if (merge_tiles is False):
            new = True
            new_pan = True

        ## start the conversion
        ## write geometry
        if (output_geometry):
            if verbosity > 1: print('Reading per pixel geometry')
            if geometry_type == 'grids': ## default s2 5x5 km grids
                xnew = np.linspace(0, grmeta['VIEW']['Average_View_Zenith'].shape[1]-1, int(global_dims[1]))
                ynew = np.linspace(0, grmeta['VIEW']['Average_View_Zenith'].shape[0]-1, int(global_dims[0]))
                if limit is not None:
                    if (dct_prj['proj4_string'] == dct['proj4_string']):
                        xnew=xnew[sub[0]:sub[0]+sub[2]]
                        ynew=ynew[sub[1]:sub[1]+sub[3]]
                    else:
                        stop
                sza = ac.shared.tiles_interp(grmeta['SUN']['Zenith'], xnew, ynew, smooth=False, method='linear')
                saa = ac.shared.tiles_interp(grmeta['SUN']['Azimuth'], xnew, ynew, smooth=False, method='linear')
                vza = ac.shared.tiles_interp(grmeta['VIEW']['Average_View_Zenith'], xnew, ynew, smooth=False, method='nearest')
                vaa = ac.shared.tiles_interp(grmeta['VIEW']['Average_View_Azimuth'], xnew, ynew, smooth=False, method='nearest')
                mask = (vaa == 0) * (vza == 0) * (saa == 0) * (sza == 0)
            elif geometry_type == 'gpt': ## use snap gpt to get nicer angles
                geometry_parameters = ['view_zenith_mean','view_azimuth_mean','sun_zenith','sun_azimuth']
                geometry_files = ac.sentinel2.gpt_geometry(bundle, output=output, target_res=geometry_res,
                                                           verbosity=verbosity, override=geometry_override)
                if geometry_format == 'GeoTIFF':
                    szai = [i for i, f in enumerate(geometry_files) if 'sun_zenith' in f][0]
                    saai = [i for i, f in enumerate(geometry_files) if 'sun_azimuth' in f][0]
                    vzai = [i for i, f in enumerate(geometry_files) if 'view_zenith_mean' in f][0]
                    vaai = [i for i, f in enumerate(geometry_files) if 'view_azimuth_mean' in f][0]
                    sza = ac.shared.read_band(geometry_files[szai], sub=sub, warp_to=warp_to)
                    saa = ac.shared.read_band(geometry_files[saai], sub=sub, warp_to=warp_to)
                    vza = ac.shared.read_band(geometry_files[vzai], sub=sub, warp_to=warp_to)
                    vaa = ac.shared.read_band(geometry_files[vaai], sub=sub, warp_to=warp_to)
                    mask = (vza == 0) * (vaa == 0)
            vza[mask] = np.nan
            sza[mask] = np.nan
            raa = (saa-vaa)
            ## negative raa
            tmp = np.where(raa<0)
            raa[tmp]=np.abs(raa[tmp])
            ## raa along 180 degree symmetry
            tmp = np.where(raa>180)
            raa[tmp]=np.abs(raa[tmp] - 360)
            raa[mask] = np.nan
            vaa = None
            saa = None
            mask = None
            ac.output.nc_write(ofile, 'raa', raa, replace_nan=True, attributes=gatts, new=new)
            if verbosity > 1: print('Wrote raa')
            new = False
            ac.output.nc_write(ofile, 'vza', vza, replace_nan=True)
            if verbosity > 1: print('Wrote vza')
            ac.output.nc_write(ofile, 'sza', sza, replace_nan=True)
            if verbosity > 1: print('Wrote sza')
            sza = None
            vza = None

        ## write lat/lon
        if (output_geolocation):
            if (os.path.exists(ofile) & (not new)):
                datasets = ac.shared.nc_datasets(ofile)
            else:
                datasets = []
            if ('lat' not in datasets) or ('lon' not in datasets):
                if verbosity > 1: print('Writing geolocation lon/lat')
                lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel=True)
                ac.output.nc_write(ofile, 'lon', lon, attributes=gatts, new=new, double=True)
                if verbosity > 1: print('Wrote lon')
                ac.output.nc_write(ofile, 'lat', lat, double=True)
                if verbosity > 1: print('Wrote lat')
                new=False

        ## write x/y
        if (output_xy):
            if os.path.exists(ofile) & (not new):
                datasets = ac.shared.nc_datasets(ofile)
            else:
                datasets = []
            if ('x' not in datasets) or ('y' not in datasets):
                if verbosity > 1: print('Writing geolocation x/y')
                x, y = ac.shared.projection_geo(dct_prj, xy=True, add_half_pixel=True)
                ac.output.nc_write(ofile, 'x', x, new=new)
                if verbosity > 1: print('Wrote x')
                ac.output.nc_write(ofile, 'y', y)
                if verbosity > 1: print('Wrote y')
                new=False

        ## write TOA bands
        quant = float(meta['QUANTIFICATION_VALUE'])
        if verbosity > 1: print('Converting bands')
        for b in rsr_bands:
            Bn = 'B{}'.format(b)
            if Bn not in safe_files[granule]: continue
            if os.path.exists(safe_files[granule][Bn]['path']):
                if b in waves_names:
                    data = ac.shared.read_band(safe_files[granule][Bn]['path'], sub=sub, warp_to=warp_to)
                    data = data.astype(np.float32)/quant
                    data[np.where(data == 0.0)] = np.nan
                    ds = 'rhot_{}'.format(waves_names[b])
                    ds_att = {'wavelength':waves_mu[b]*1000}
                    #for k in band_data: ds_att[k] = band_data[k][b]
                    if percentiles_compute:
                        ds_att['percentiles'] = percentiles
                        ds_att['percentiles_data'] = np.nanpercentile(data, percentiles)
                    ## write to ms file
                    ac.output.nc_write(ofile, ds, data, replace_nan=True, attributes=gatts, new=new, dataset_attributes = ds_att)
                    new = False
                    if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data.shape))
            else:
                continue

        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))

        if limit is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles)
