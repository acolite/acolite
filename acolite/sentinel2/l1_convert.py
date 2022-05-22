## def l1_convert
## converts sentinel .SAFE bundle data to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-02-11
## modifications: 2021-10-13 (QV) support for new L1C format from processing baseline 4
#                 2021-10-14 (QV) fixed band specific footprints for band specific geometry for PB004
##                2021-12-08 (QV) added nc_projection
##                2021-12-31 (QV) new handling of settings

def l1_convert(inputfile, output = None, settings = {},
                percentiles_compute = True,
                percentiles = (0,1,5,10,25,50,75,90,95,99,100),
                check_sensor = True,
                check_time = True,
                max_merge_time = 600, # seconds
                geometry_format='GeoTIFF', ## for gpt geometry
                verbosity = 5):

    import sys, os, glob, dateutil.parser, time
    from osgeo import ogr,osr,gdal
    import acolite as ac
    import scipy.ndimage
    import numpy as np
    t0 = time.time()

    if 'verbosity' in settings: verbosity = settings['verbosity']

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    new = True
    warp_to = None

    ofile = None
    ofiles = []
    for bundle in inputfile:
        #if output is None: output = os.path.dirname(bundle)
        if verbosity > 1: print('Starting conversion of {}'.format(bundle))

        try:
            safe_files = ac.sentinel2.safe_test(bundle)
        except:
            print('File not recognised: {}'.format(bundle))
            #print("Error:", sys.exc_info()[0])
            continue

        if 'granules' not in safe_files:
            print('File not recognised: {}'.format(bundle))
            continue

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
            print('Processing of {} Sentinel-2 {} data not supported'.format(bundle, meta['PROCESSING_LEVEL']))
            continue

        ## merge sensor specific settings
        if new:
            setu = ac.acolite.settings.parse(sensor, settings=settings)
            verbosity = setu['verbosity']
            if output is None: output = setu['output']

            gains = setu['gains']
            gains_toa = setu['gains_toa']
            offsets = setu['offsets']
            offsets_toa = setu['offsets_toa']

            s2_target_res=setu['s2_target_res']

            geometry_type=setu['geometry_type']
            geometry_res=setu['geometry_res']
            geometry_per_band=setu['geometry_per_band']
            geometry_fixed_footprint=setu['geometry_fixed_footprint']

            s2_include_auxillary = setu['s2_include_auxillary']
            s2_project_auxillary = setu['s2_project_auxillary']
            netcdf_projection = setu['netcdf_projection']

            output_geolocation=setu['output_geolocation']
            output_xy=setu['output_xy']
            output_geometry=setu['output_geometry']
            vname = setu['region_name']

            limit=setu['limit']
            ## check if ROI polygon is given
            if setu['polylakes']:
                poly = ac.shared.polylakes(setu['polylakes_database'])
                setu['polygon_limit'] = False
            else:
                poly = setu['polygon']
            clip, clip_mask = False, None
            if poly is not None:
                if os.path.exists(poly):
                    try:
                        limit = ac.shared.polygon_limit(poly)
                        if setu['polygon_limit']:
                            print('Using limit from polygon envelope: {}'.format(limit))
                        else:
                            limit = setu['limit']
                        clip = True
                    except:
                        print('Failed to import polygon {}'.format(poly))

            ## check if merging settings make sense
            merge_tiles = setu['merge_tiles']
            merge_zones = setu['merge_zones']
            extend_region = setu['extend_region']
            if (limit is None) & (merge_tiles):
                if verbosity > 0: print("Merging tiles not supported without ROI limit")
                merge_tiles = False
            if merge_tiles:
                merge_zones = True
                extend_region = True
        sub = None

        dtime = dateutil.parser.parse(grmeta['SENSING_TIME'])
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()
        dct = ac.sentinel2.projection(grmeta, s2_target_res=int(s2_target_res))
        global_dims = dct['dimensions']

        #mgrs_tile = meta['PRODUCT_URI'].split('_')[-2]
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

        ## gains
        gains_dict = None
        if gains & (gains_toa is not None):
            if len(gains_toa) == len(rsr_bands):
                gains_dict = {b: float(gains_toa[ib]) for ib, b in enumerate(rsr_bands)}

        ## offsets
        offsets_dict = None
        if offsets & (offsets_toa is not None):
            if len(offsets_toa) == len(rsr_bands):
                offsets_dict = {b: float(offsets_toa[ib]) for ib, b in enumerate(rsr_bands)}

        ## get F0 - not stricty necessary if using USGS reflectance
        f0 = ac.shared.f0_get()
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data'])*10, rsr)

        ## make global attributes for L1R NetCDF
        gatts = {'sensor':sensor, 'isodate':isodate, 'global_dims':global_dims,
                 'sza':sza, 'vza':vza, 'raa':raa, 'vaa': vaa, 'saa': saa, 'se_distance': se_distance,
                 'mus': np.cos(sza*(np.pi/180.)), 'granule': granule, 'mgrs_tile': mgrs_tile,
                 'acolite_file_type': 'L1R'}
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

        ## get projection info for netcdf
        if netcdf_projection:
            nc_projection = ac.shared.projection_netcdf(dct_prj, add_half_pixel=True)
        else:
            nc_projection = None

        ## save projection keys in gatts
        pkeys = ['xrange', 'yrange', 'proj4_string', 'pixel_size', 'zone']
        for k in pkeys:
            if k in dct_prj: gatts[k] = dct_prj[k]

        ## with subsetting fix the offsets should not be required 2021-10-28
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

        ## if we are clipping to a given polygon get the clip_mask here
        if clip:
            clip_mask = ac.shared.polygon_crop(dct_prj, poly, return_sub=False)
            clip_mask = clip_mask.astype(bool) == False
            print('clip mask', clip_mask.shape)

        ## start the conversion
        ## write geometry
        if (output_geometry):
            if verbosity > 1: print('Reading per pixel geometry')
            if (geometry_type == 'grids') | (geometry_type == 'grids_footprint'):
                ## just use 60m band for geometry, it will be interpolated later
                if geometry_res == 10:
                    target_file = safe_files[granule]['B2']['path']
                elif geometry_res == 20:
                    target_file = safe_files[granule]['B11']['path']
                elif geometry_res == 60:
                    target_file = safe_files[granule]['B1']['path']

                ## get proj dct for geometry
                dct_geom = ac.sentinel2.projection(grmeta, s2_target_res=geometry_res)
                ## with subsetting fix the offsets should not be required 2021-10-28
                xyr_geom = [min(dct_geom['xrange']),
                            min(dct_geom['yrange']),
                            max(dct_geom['xrange']),
                            max(dct_geom['yrange']),
                            dct_geom['proj4_string']]

                ## warp settings for read_band
                warp_to_geom = (dct_geom['proj4_string'], xyr_geom, dct_geom['pixel_size'][0], dct_geom['pixel_size'][1], 'average')

                ## open target file to get dimensions
                g = gdal.Open(target_file)
                xSrc = g.RasterXSize
                ySrc = g.RasterYSize
                g = None
                xnew = np.linspace(0, grmeta['VIEW']['Average_View_Zenith'].shape[1]-1, int(xSrc))
                ynew = np.linspace(0, grmeta['VIEW']['Average_View_Zenith'].shape[0]-1, int(ySrc))
                sza = ac.shared.tiles_interp(grmeta['SUN']['Zenith'], xnew, ynew, smooth=False, method='linear')
                saa = ac.shared.tiles_interp(grmeta['SUN']['Azimuth'], xnew, ynew, smooth=False, method='linear')

                ## default s2 5x5 km grids
                if geometry_type == 'grids':
                    #xnew = np.linspace(0, grmeta['VIEW']['Average_View_Zenith'].shape[1]-1, int(global_dims[1]))
                    #ynew = np.linspace(0, grmeta['VIEW']['Average_View_Zenith'].shape[0]-1, int(global_dims[0]))
                    vza = ac.shared.tiles_interp(grmeta['VIEW']['Average_View_Zenith'], xnew, ynew, smooth=False, method='nearest')
                    vaa = ac.shared.tiles_interp(grmeta['VIEW']['Average_View_Azimuth'], xnew, ynew, smooth=False, method='nearest')

                ## use s2 5x5 km grids with detector footprint interpolation
                if geometry_type == 'grids_footprint':
                    ## compute vza and saa
                    gml_files = glob.glob('{}/GRANULE/{}/QI_DATA/*MSK_DETFOO*.gml'.format(bundle, granule))
                    gml_files.sort()

                    jp2_files = glob.glob('{}/GRANULE/{}/QI_DATA/*MSK_DETFOO*.jp2'.format(bundle, granule))
                    jp2_files.sort()

                    ## get detector footprint for 10/20/60 m band
                    if len(gml_files) > 0:
                        dval, dfoo = ac.sentinel2.detector_footprint(target_file, gml_files[0])
                    elif len(jp2_files) > 0:
                        dfoo = ac.shared.read_band(jp2_files[0], warp_to=warp_to_geom)
                        dval = np.unique(dfoo)
                    else:
                        print('No footprint files found')
                        continue
                    bands = [str(bi) for bi, b in enumerate(rsr_bands)]

                    ## set vza and vaa size to geometry res target size
                    vza = np.zeros((int(dfoo.shape[0]), int(dfoo.shape[1])))+np.nan
                    vaa = np.zeros((int(dfoo.shape[0]), int(dfoo.shape[1])))+np.nan

                    if verbosity>1:print('Computing band average per detector geometry')
                    for nf, bv in enumerate(dval):
                        if bv == 0: continue ## fill value in new format
                        ## compute detector average geometry
                        if verbosity>2:print('Detector {}'.format(bv))
                        ave_vza = None
                        ave_vaa = None
                        for b in bands:
                            if '{}'.format(bv) not in grmeta['VIEW_DET'][b]: continue
                            bza = grmeta['VIEW_DET'][b]['{}'.format(bv)]['Zenith']
                            baa = grmeta['VIEW_DET'][b]['{}'.format(bv)]['Azimuth']
                            bza = ac.sentinel2.grid_extend(bza, iterations=1, crop=False)
                            baa = ac.sentinel2.grid_extend(baa, iterations=1, crop=False)
                            if ave_vaa is None:
                                ave_vza = bza
                                ave_vaa = baa
                            else:
                                ave_vza = np.dstack((ave_vza, bza))
                                ave_vaa = np.dstack((ave_vaa, baa))

                        if ave_vza is not None:
                            ave_vza = np.nanmean(ave_vza, axis=2)
                            ave_vaa = np.nanmean(ave_vaa, axis=2)
                            ## end compute detector average geometry
                            ## interpolate grids to current detector
                            det_mask = dfoo==bv
                            ## add +1 to xnew and ynew since we are not cropping the extended grid
                            vza[det_mask] = ac.shared.tiles_interp(ave_vza, xnew+1, ynew+1, smooth=False, fill_nan=True,
                                                          target_mask = det_mask, target_mask_full=False, method='linear')
                            vaa[det_mask] = ac.shared.tiles_interp(ave_vaa, xnew+1, ynew+1, smooth=False, fill_nan=True,
                                                          target_mask = det_mask, target_mask_full=False, method='linear')

                ## use target band so we can just do the 60 metres geometry
                if os.path.exists(target_file):
                    sza = ac.shared.warp_from_source(target_file, dct_prj, sza, warp_to=warp_to) # alt (dct, dct_prj, sza)
                    saa = ac.shared.warp_from_source(target_file, dct_prj, saa, warp_to=warp_to)
                    vza = ac.shared.warp_from_source(target_file, dct_prj, vza, warp_to=warp_to)
                    vaa = ac.shared.warp_from_source(target_file, dct_prj, vaa, warp_to=warp_to)
                    mask = (vaa == 0) * (vza == 0) * (saa == 0) * (sza == 0)
                else:
                    print('Could not access {}'.format(target_file))
                    print('Path length {} greater than the recommended path length on Windows'.format(len(target_file)))
                    return(ofiles, setu)

                ## compute band specific geometry
                if geometry_per_band:
                    print('Computing band specific per detector geometry')
                    vza_all = np.zeros((int(dfoo.shape[0]), int(dfoo.shape[1]), len(bands)))+np.nan
                    vaa_all = np.zeros((int(dfoo.shape[0]), int(dfoo.shape[1]), len(bands)))+np.nan

                    grid_shape = grmeta['VIEW']['0']['Zenith'].shape[0]+2, grmeta['VIEW']['0']['Zenith'].shape[1]+2
                    vza_grid = np.zeros((grid_shape[0], grid_shape[1], len(bands)))+np.nan
                    vaa_grid = np.zeros((grid_shape[0], grid_shape[1], len(bands)))+np.nan

                    ## use footprint from B1
                    if geometry_fixed_footprint:
                        gml_files = glob.glob('{}/GRANULE/{}/QI_DATA/*MSK_DETFOO*.gml'.format(bundle, granule))
                        gml_files.sort()
                        jp2_files = glob.glob('{}/GRANULE/{}/QI_DATA/*MSK_DETFOO*.jp2'.format(bundle, granule))
                        jp2_files.sort()

                        ## get detector footprint for 10/20/60 m band
                        if len(gml_files) > 0:
                            dval, dfoo = ac.sentinel2.detector_footprint(target_file, gml_files[0])
                        elif len(jp2_files) > 0:
                            dfoo = ac.shared.read_band(jp2_files[0], warp_to=warp_to_geom)
                            dval = np.unique(dfoo)

                    ## compute band specific view geometry
                    for bi, b in enumerate(bands):
                        Bn = band_data['BandNames'][b]
                        print('Computing band specific geometry - {}'.format(Bn))

                        ## band specific footprint
                        if not geometry_fixed_footprint:
                            gml = glob.glob('{}/GRANULE/{}/QI_DATA/*MSK_DETFOO_B{}.gml'.format(bundle, granule, Bn[1:].zfill(2)))
                            jp2 = glob.glob('{}/GRANULE/{}/QI_DATA/*MSK_DETFOO_B{}.jp2'.format(bundle, granule, Bn[1:].zfill(2)))

                            if len(gml) > 0:
                                dval, dfoo = ac.sentinel2.detector_footprint(target_file, gml[0])
                            elif len(jp2) > 0:
                                dfoo = ac.shared.read_band(jp2[0], warp_to=warp_to_geom)
                                dval = np.unique(dfoo)

                        for nf, bv in enumerate(dval):
                            if bv == 0: continue ## fill value in new format
                            print('Computing band specific geometry - {} Detector {}'.format(Bn, bv))
                            det_mask = dfoo==bv

                            bza = grmeta['VIEW_DET'][b]['{}'.format(bv)]['Zenith']
                            baa = grmeta['VIEW_DET'][b]['{}'.format(bv)]['Azimuth']
                            bza = ac.sentinel2.grid_extend(bza, iterations=1, crop=False)
                            baa = ac.sentinel2.grid_extend(baa, iterations=1, crop=False)

                            ## add detector to band VZA grid
                            vza_tmp =  vza_grid[:,:,bi]
                            ang_sub = np.where(np.isfinite(bza))
                            vza_tmp[ang_sub] = bza[ang_sub]
                            vza_grid[:,:,bi] = vza_tmp

                            ## add detector to band VAA grid
                            vaa_tmp =  vaa_grid[:,:,bi]
                            ang_sub = np.where(np.isfinite(baa))
                            vaa_tmp[ang_sub] = baa[ang_sub]
                            vaa_grid[:,:,bi] = vaa_tmp

                            ## add +1 to xnew and ynew since we are not cropping the extended grid
                            vza_tmp = vza_all[:,:, bi]
                            vza_tmp[det_mask] = ac.shared.tiles_interp(vza_grid[:,:,bi], xnew+1, ynew+1, smooth=False, fill_nan=True,
                                                                          target_mask = det_mask, target_mask_full=False, method='linear')
                            vza_all[:,:, bi] = vza_tmp

                            vaa_tmp = vaa_all[:,:, bi]
                            vaa_tmp[det_mask] = ac.shared.tiles_interp(vaa_grid[:,:,bi], xnew+1, ynew+1, smooth=False, fill_nan=True,
                                                                          target_mask = det_mask, target_mask_full=False, method='linear')
                            vaa_all[:,:, bi] = vaa_tmp

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

            ## add out of swath masks to vza and sza
            vza[mask] = np.nan
            sza[mask] = np.nan

            ## clip geometry data
            if clip:
                sza[clip_mask] = np.nan
                saa[clip_mask] = np.nan
                vza[clip_mask] = np.nan
                vaa[clip_mask] = np.nan

            if setu['s2_write_vaa']:
                ac.output.nc_write(ofile, 'vaa', vaa, replace_nan=True,
                                        attributes=gatts, new=new, nc_projection=nc_projection,
                                                netcdf_compression=setu['netcdf_compression'],
                                                netcdf_compression_level=setu['netcdf_compression_level'])
                new = False
                if verbosity > 1: print('Wrote vaa {}'.format(vaa.shape))

            if setu['s2_write_saa']:
                ac.output.nc_write(ofile, 'saa', saa, replace_nan=True,
                                        attributes=gatts, new=new, nc_projection=nc_projection,
                                                netcdf_compression=setu['netcdf_compression'],
                                                netcdf_compression_level=setu['netcdf_compression_level'])
                new = False
                if verbosity > 1: print('Wrote saa {}'.format(saa.shape))

            ## compute relative azimuth angle
            raa = np.abs(saa-vaa)
            ## raa along 180 degree symmetry
            tmp = np.where(raa>180)
            raa[tmp]=np.abs(raa[tmp] - 360)
            raa[mask] = np.nan
            vaa = None

            ac.output.nc_write(ofile, 'raa', raa, replace_nan=True,
                                    attributes=gatts, new=new, nc_projection=nc_projection,
                                                netcdf_compression=setu['netcdf_compression'],
                                                netcdf_compression_level=setu['netcdf_compression_level'])
            if verbosity > 1: print('Wrote raa {}'.format(raa.shape))
            raa = None
            new = False
            ac.output.nc_write(ofile, 'vza', vza, replace_nan=True,
                                netcdf_compression=setu['netcdf_compression'],
                                netcdf_compression_level=setu['netcdf_compression_level'])
            if verbosity > 1: print('Wrote vza {}'.format(vza.shape))
            ac.output.nc_write(ofile, 'sza', sza, replace_nan=True,
                                netcdf_compression=setu['netcdf_compression'],
                                netcdf_compression_level=setu['netcdf_compression_level'])
            if verbosity > 1: print('Wrote sza {}'.format(sza.shape))
            sza = None
            vza = None

            ## write per band geometry
            if (geometry_per_band) & ((geometry_type == 'grids') | (geometry_type == 'grids_footprint')):
                for bi, b in enumerate(rsr_bands):
                    Bn = 'B{}'.format(b)
                    print('Writing view geometry for {} {} nm'.format(Bn, waves_names[b]))
                    ## band specific view zenith angle
                    vza = ac.shared.warp_from_source(target_file, dct_prj, vza_all[:,:,bi], warp_to=warp_to)
                    ac.output.nc_write(ofile, 'vza_{}'.format(waves_names[b]), vza, replace_nan=True,
                                        netcdf_compression=setu['netcdf_compression'],
                                        netcdf_compression_level=setu['netcdf_compression_level'])
                    vza = None
                    ## band specific view azimuth angle
                    vaa = ac.shared.warp_from_source(target_file, dct_prj, vaa_all[:,:,bi], warp_to=warp_to)
                    if setu['s2_write_vaa']:
                        ac.output.nc_write(ofile, 'vaa_{}'.format(waves_names[b]), vaa, replace_nan=True,
                                                attributes=gatts, new=new, nc_projection=nc_projection,
                                                        netcdf_compression=setu['netcdf_compression'],
                                                        netcdf_compression_level=setu['netcdf_compression_level'])
                    ## compute relative azimuth angle
                    raa = np.abs(saa-vaa)
                    vaa = None
                    tmp = np.where(raa>180)
                    raa[tmp]=np.abs(raa[tmp] - 360)
                    raa[mask] = np.nan
                    ac.output.nc_write(ofile, 'raa_{}'.format(waves_names[b]), raa, replace_nan=True,
                                       netcdf_compression=setu['netcdf_compression'],
                                       netcdf_compression_level=setu['netcdf_compression_level'])
                    raa = None
            ## delete sun azimuth & mask
            saa = None
            mask = None

        ## write lat/lon
        if (output_geolocation):
            if (os.path.exists(ofile) & (not new)):
                datasets = ac.shared.nc_datasets(ofile)
            else:
                datasets = []
            if ('lat' not in datasets) or ('lon' not in datasets):
                if verbosity > 1: print('Writing geolocation lon/lat')
                lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel=True)
                print(lon.shape)
                ac.output.nc_write(ofile, 'lon', lon, attributes=gatts, new=new, double=True,
                                    nc_projection=nc_projection,
                                    netcdf_compression=setu['netcdf_compression'],
                                    netcdf_compression_level=setu['netcdf_compression_level'])
                lon = None
                if verbosity > 1: print('Wrote lon')
                print(lat.shape)
                ac.output.nc_write(ofile, 'lat', lat, double=True,
                                    netcdf_compression=setu['netcdf_compression'],
                                    netcdf_compression_level=setu['netcdf_compression_level'])
                lat = None
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
                ac.output.nc_write(ofile, 'x', x, new=new,
                                    netcdf_compression=setu['netcdf_compression'],
                                    netcdf_compression_level=setu['netcdf_compression_level'],
                                    netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
                x = None
                if verbosity > 1: print('Wrote x')
                ac.output.nc_write(ofile, 'y', y,
                                    netcdf_compression=setu['netcdf_compression'],
                                    netcdf_compression_level=setu['netcdf_compression_level'],
                                    netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
                y = None
                if verbosity > 1: print('Wrote y')
                new=False

        ## auxillary data
        if s2_include_auxillary:
            ofile_aux = '{}/{}'.format(os.path.dirname(ofile), os.path.basename(ofile).replace('_L1R.nc', '_AUX.nc'))
            ofile_aux_new = True
            for source in ['AUX_CAMSFO', 'AUX_ECMWFT']:
                ## read aux data
                aux_data = ac.sentinel2.auxillary(bundle, granule, sources=[source])
                if len(aux_data) > 0:
                    ## add to gatts
                    for ai, an in enumerate(aux_data):
                        gatts['{}_{}_{}'.format(source, an, 'dimensions')] = aux_data[an]['values'].shape
                        gatts['{}_{}_{}'.format(source, an, 'values')] = aux_data[an]['values'].flatten()
                        gatts['{}_{}_{}'.format(source, an, 'longitudes')] = aux_data[an]['longitudes'].flatten()
                        gatts['{}_{}_{}'.format(source, an, 'latitudes')] = aux_data[an]['latitudes'].flatten()
                    ## project to extent
                    if s2_project_auxillary:
                        aux_file = '{}/GRANULE/{}/AUX_DATA/{}'.format(bundle, granule, source)
                        # gdal warp
                        #adata = ac.shared.read_band(aux_file, sub=None, warp_to=warp_to)
                        lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel=True)
                        llo = np.vstack((lon.flatten(),lat.flatten())).T
                        lon = None
                        lat = None
                        for ai, an in enumerate(aux_data):
                            lli = np.stack((aux_data[an]['longitudes'].flatten(),
                                            aux_data[an]['latitudes'].flatten())).T
                            v = aux_data[an]['values'].flatten()
                            ## interpolate and fill edges
                            ret = scipy.interpolate.griddata(lli, v, llo)
                            ret = ac.shared.fillnan(ret.reshape(int(gatts['global_dims'][0]), int(gatts['global_dims'][1])))
                            ## write
                            ac.output.nc_write(ofile_aux, '{}_{}'.format(source, an), ret, replace_nan=True,
                                                attributes=gatts, new = ofile_aux_new, nc_projection=nc_projection,
                                                netcdf_compression=setu['netcdf_compression'],
                                                netcdf_compression_level=setu['netcdf_compression_level'],
                                                netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
                            if verbosity > 1: print('Wrote {}'.format('{}_{}'.format(source, an)))
                            ret = None
                            ofile_aux_new = False

        ## write TOA bands
        quant = float(meta['QUANTIFICATION_VALUE'])
        nodata = int(meta['NODATA'])
        ## new offset in processing baseline 4
        if 'RADIO_ADD_OFFSET' in band_data:
            for Bn in band_data['RADIO_ADD_OFFSET']:
                band_data['RADIO_ADD_OFFSET'][Bn] = float(band_data['RADIO_ADD_OFFSET'][Bn])
        if verbosity > 1: print('Converting bands')
        for bi, b in enumerate(rsr_bands):
            Bn = 'B{}'.format(b)
            if Bn not in safe_files[granule]: continue
            if os.path.exists(safe_files[granule][Bn]['path']):
                if b in waves_names:
                    data = ac.shared.read_band(safe_files[granule][Bn]['path'], sub=sub, warp_to=warp_to)
                    data_mask = data == nodata
                    data = data.astype(np.float32)
                    if 'RADIO_ADD_OFFSET' in band_data: data += band_data['RADIO_ADD_OFFSET'][Bn]
                    data /= quant
                    data[data_mask] = np.nan
                    if clip: data[clip_mask] = np.nan
                    ds = 'rhot_{}'.format(waves_names[b])
                    ds_att = {'wavelength':waves_mu[b]*1000}

                    if gains & (gains_dict is not None):
                        ds_att['toa_gain'] = gains_dict[b]
                        data *= ds_att['toa_gain']
                        if verbosity > 1: print('Converting bands: Applied TOA gain {} to {}'.format(ds_att['toa_gain'], ds))
                    if offsets & (offsets_dict is not None):
                        ds_att['toa_offset'] = offsets_dict[b]
                        data *= ds_att['toa_offset']
                        if verbosity > 1: print('Converting bands: Applied TOA offset {} to {}'.format(ds_att['toa_gain'], ds))

                    #for k in band_data: ds_att[k] = band_data[k][b]
                    if percentiles_compute:
                        ds_att['percentiles'] = percentiles
                        ds_att['percentiles_data'] = np.nanpercentile(data, percentiles)

                    ## write to ms file
                    ac.output.nc_write(ofile, ds, data, replace_nan=True, attributes=gatts, new=new,
                                        dataset_attributes = ds_att, nc_projection=nc_projection,
                                        netcdf_compression=setu['netcdf_compression'],
                                        netcdf_compression_level=setu['netcdf_compression_level'],
                                        netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
                    new = False
                    if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data.shape))
            else:
                continue

        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))

        if limit is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
