## def l1_convert
## converts Sentinel-2 zarr to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2025-11-06
## modifications: 2025-11-17 (QV) added tile merging
##                2026-01-14 (QV) changed AUX interpolation and integrated to sentinel2.zarr
##                2026-02-18 (QV) added plane_fit_geom as keyword - but don't use it!

def l1_convert(inputfile, output = None, settings = None,
                check_sensor = True,
                check_time = True,
                max_merge_time = 600, # seconds,
                plane_fit_geom = False,
                ):

    import sys, os, glob, dateutil.parser, time
    import pyproj, zarr

    from osgeo import ogr,osr,gdal
    import acolite as ac
    import scipy.ndimage
    import numpy as np
    t0 = time.time()

    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    verbosity = setu['verbosity']

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
    ofile_aux_new = True

    ofile = None
    ofiles = []
    for bundle in inputfile:
        #if output is None: output = os.path.dirname(bundle)
        if verbosity > 1: print('Starting conversion of {}'.format(bundle))

        t0 = time.time() ## track time

        ## track if we are dealing with local or remote dataset
        if os.path.exists(bundle):
            url = False
        else:
            url = True

        ## open zarr bundle
        z = zarr.open(bundle, mode='r')

        ## get zarr metadata
        meta = ac.zarr.meta_parse(z)

        ## sensor settings
        sensor = meta['sensor']

        if meta['processing:level'] != 'L1C':
            print('Processing of {} Sentinel-2 {} data not supported'.format(bundle, meta['processing:level']))
            continue

        ## merge sensor specific settings
        if new:
            ## get sensor specific defaults
            setd = ac.acolite.settings.parse(sensor)
            ## set sensor default if user has not specified the setting
            for k in setd:
                if k not in ac.settings['user']: setu[k] = setd[k]
            ## end set sensor specific defaults

            verbosity = setu['verbosity']
            if output is None: output = setu['output']

            ## check if merging settings make sense
            extend_region = setu['extend_region']
            if setu['merge_tiles']:
                if (setu['limit'] is None):
                    if not setu['merge_full_tiles']:
                        if verbosity > 0: print("Merging tiles without ROI limit, merging to first tile extent")
                    else:
                        if verbosity > 0: print("Merging tiles without ROI limit, merging to all tiles extent")
                        dct_tiles = ac.sentinel2.zarr.multi_tile_extent(inputfile, dct = None, s2_target_res = setu['s2_target_res'])
                else:
                    extend_region = True

        ## sub is set to None
        sub = None

        ## datetime
        dtime = dateutil.parser.parse(meta['start_datetime'])
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()
        #se_distance = meta['reflectance_correction_factor_from_the_Sun_Earth_distance_variation_computed_using_the_acquisition_date']

        ## get MGRS tile info
        bn, ext = os.path.splitext(os.path.basename(bundle))
        sp = bn.split('_')
        mgrs_tile = sp[-2]

        ## get scene projection and extent
        dct = ac.sentinel2.zarr.projection(z, meta = meta, s2_target_res = setu['s2_target_res'])
        source_srs = dct['projection']

        ## grid metadata
        angle_order = z['conditions']['geometry']['angle'][:]
        detectors = z['conditions']['geometry']['detector'][:]
        bands = z['conditions']['geometry']['band'][:]

        ## angle grids
        sun_angles = z['conditions']['geometry']['sun_angles'][:]
        viewing_incidence_angles = z['conditions']['geometry']['viewing_incidence_angles'][:]

        ## mean angle grids
        mean_sun_angles = z['conditions']['geometry']['mean_sun_angles'][:]
        mean_viewing_incidence_angles = z['conditions']['geometry']['mean_viewing_incidence_angles'][:]

        ## grid x and y spacing
        x_grid = z['conditions']['geometry']['x'][:]
        y_grid = z['conditions']['geometry']['y'][:]
        ## set up grid mesh
        x_grid_mesh, y_grid_mesh = np.meshgrid(x_grid, y_grid, indexing='ij')

        ## zenith and azimuth indices
        zi = np.where(angle_order == 'zenith')[0][0]
        ai = np.where(angle_order == 'azimuth')[0][0]

        ## get mean geometry
        sza = mean_sun_angles[zi]
        saa = mean_sun_angles[ai]
        vza = np.nanmean(mean_viewing_incidence_angles[:, zi])
        vaa = np.nanmean(mean_viewing_incidence_angles[:, ai])
        raa = np.abs(saa-vaa)
        while raa > 180: raa = np.abs(360 - raa)

        ## read RSR
        if setu['rsr_version'] is not None:
            sensor_lut = '{}_{}'.format(sensor, setu['rsr_version'])
        else:
            sensor_lut = '{}'.format(sensor)
        rsrd = ac.shared.rsr_dict(sensor_lut)[sensor_lut]
        waves_names = rsrd['wave_mu']
        rsr_bands = rsrd['rsr_bands']

        ## parse gains
        gains_dict = None
        if setu['gains'] & (setu['gains_toa'] is not None):
            if len(setu['gains_toa']) == len(rsr_bands):
                gains_dict = {b: float(setu['gains_toa'][ib]) for ib, b in enumerate(rsr_bands)}

        ## parse offsets
        offsets_dict = None
        if setu['offsets'] & (setu['offsets_toa'] is not None):
            if len(setu['offsets_toa']) == len(rsr_bands):
                offsets_dict = {b: float(setu['offsets_toa'][ib]) for ib, b in enumerate(rsr_bands)}

        ## get F0
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data'])*10, rsrd['rsr'])

        ## make global attributes for L1R NetCDF
        gatts = {'sensor':sensor, 'isodate':isodate, 'global_dims': dct['dimensions'],
                 'sza':sza, 'vza':vza, 'raa':raa, 'vaa': vaa, 'saa': saa, 'se_distance': se_distance,
                 'mus': np.cos(sza*(np.pi/180.)), 'granule': bn, 'mgrs_tile': mgrs_tile,
                 'acolite_file_type': 'L1R'}
        if setu['merge_tiles']:
            gatts['tile_code'] = 'merged'
        else:
            gatts['tile_code'] = '{}'.format(gatts['mgrs_tile'])

        stime = dateutil.parser.parse(gatts['isodate'])
        oname = '{}_{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'), gatts['tile_code'])
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ## output file information
        if (setu['merge_tiles'] is False) | (ofile is None):
            ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
            gatts['oname'] = oname
            gatts['ofile'] = ofile
        elif (setu['merge_tiles']) & (ofile is None):
            ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
            gatts['oname'] = oname
            gatts['ofile'] = ofile

        ## check if we should merge these tiles
        if (setu['merge_tiles']) & (not new) & (os.path.exists(ofile)):
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
            gatts['{}_wave'.format(b)] = rsrd['wave_nm'][b]
            gatts['{}_name'.format(b)] = rsrd['wave_name'][b]
            gatts['{}_f0'.format(b)] = f0_b[b]

        ## full scene
        gatts['scene_xrange'] = dct['xrange']
        gatts['scene_yrange'] = dct['yrange']
        #gatts['scene_proj4_string'] = dct['proj4_string']
        gatts['scene_pixel_size'] = dct['pixel_size']
        gatts['scene_dims'] = dct['dimensions']
        #if 'zone' in dct: gatts['scene_zone'] = dct['zone']

        ## check crop
        if (sub is None) & (setu['limit'] is not None):
            dct_sub = ac.shared.projection_sub(dct, setu['limit'], four_corners = True)
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

        ## remove warp_to from previous run if merge_tiles is not set
        if (setu['merge_tiles'] is False): warp_to = None
        if sub is None: ## full tile processing
            ## determine warping target
            if (warp_to is None):
                if (setu['merge_tiles'] & setu['merge_full_tiles']): ## warp to all tile extent
                    dct_prj = {k:dct_tiles[k] for k in dct_tiles}
                else: ## warp to current/first tile
                    dct_prj = {k:dct[k] for k in dct}
        else:
            gatts['sub'] = sub
            gatts['limit'] = setu['limit']

            ## get the target NetCDF dimensions and dataset offset
            if (warp_to is None):
                if (extend_region): ## include part of the roi not covered by the scene
                    dct_prj = {k:dct_sub['region'][k] for k in dct_sub['region']}
                else: ## just include roi that is covered by the scene
                    dct_prj = {k:dct_sub[k] for k in dct_sub}
        ## end cropped

        ## get projection info for netcdf
        if setu['netcdf_projection']:
            nc_projection = ac.shared.projection_netcdf(dct_prj, add_half_pixel = True)
        else:
            nc_projection = None

        ## save projection keys in gatts
        pkeys = ['xrange', 'yrange', 'projection', 'pixel_size', 'zone']
        for k in pkeys:
            if k in dct_prj: gatts[k] = dct_prj[k]

        ## warp settings for read_band, using average resampling method
        warp_to = ac.shared.projection_warp_to(dct_prj, res_method = 'average')

        ## store scene and output dimensions
        gatts['scene_dims'] = dct['ydim'], dct['xdim']
        gatts['global_dims'] = dct_prj['dimensions']

        ## new file for every bundle if not merging
        if (setu['merge_tiles'] is False):
            new = True
            ofile_aux_new = True

        ## if we are clipping to a given polygon get the clip_mask here
        if setu['polygon_clip']:
            clip_mask = ac.shared.polygon_crop(dct_prj, setu['polygon'], return_sub = False)
            clip_mask = clip_mask.astype(bool) == False
            print('clip mask', clip_mask.shape)

        if new:
            gemo = ac.gem.gem(ofile, new = True)
            gemo.gatts = {k: gatts[k] for k in gatts}
            gemo.nc_projection = nc_projection
            new = False
        else:
            gatts =  {k: gemo.gatts[k] for k in gemo.gatts} ## read gatts to be updated

        ## set up output file
        if os.path.exists(ofile) & (not new):
            gemo.datasets_read()
            datasets = gemo.datasets
        else:
            datasets = []

        ## write lat/lon
        if (setu['output_geolocation']):
            if ('lat' not in datasets) or ('lon' not in datasets):
                if verbosity > 1: print('Writing geolocation lon/lat')
                lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel = True)
                gemo.write('lon', lon)
                if verbosity > 1: print('Wrote lon {}'.format(lon.shape))
                lon = None
                gemo.write('lat', lat)
                if verbosity > 1: print('Wrote lat {}'.format(lat.shape))
                lat = None

        ## write x/y
        if (setu['output_xy']):
            if ('xm' not in datasets) or ('ym' not in datasets):
                if verbosity > 1: print('Writing geolocation x/y')
                x, y = ac.shared.projection_geo(dct_prj, xy = True, add_half_pixel = True)
                gemo.write('xm', x)
                if verbosity > 1: print('Wrote xm {}'.format(x.shape))
                x = None
                gemo.write('ym', y)
                if verbosity > 1: print('Wrote ym {}'.format(y.shape))
                y = None

        # dct_prj tracks projection for target scene
        #band_x_mesh, band_y_mesh = ac.shared.projection_geo(dct_prj, xy = True, add_half_pixel = True)

        ## use dct_sub for current tile grid mesh
        #band_x_mesh, band_y_mesh = ac.shared.projection_geo(dct_sub, xy = True, add_half_pixel = True)
        #band_x_mesh, band_y_mesh = ac.shared.projection_geo(dct_sub['region'], xy = True, add_half_pixel = True)

        ## use full tile grid mesh, reproject and subset later
        ## needed for tile merging in different zones
        print('Constructing interpolator mesh')
        band_x_mesh, band_y_mesh = ac.shared.projection_geo(dct, xy = True, add_half_pixel = True)
        print('Interpolator mesh shape', band_x_mesh.shape)

        ## 60 metre warp to for geometry
        warp_to_geom = ac.shared.projection_warp_to(dct, res_method = 'average')

        ## coordinates for geometry interpolator
        xnew = np.linspace(0, x_grid.shape[0]-1, num = band_x_mesh.shape[0])
        ynew = np.linspace(0, y_grid.shape[0]-1, num = band_y_mesh.shape[1])

        ## for per pixel geometry
        if (setu['output_geometry']):
            print('Computing per pixel geometries')

            ## dct_prj tracks projection for target scene
            #band_x_mesh, band_y_mesh = ac.shared.projection_geo(dct_prj, xy = True, add_half_pixel = True)

            ## use dct_sub for current tile
            #band_x_mesh, band_y_mesh = ac.shared.projection_geo(dct_sub, xy = True, add_half_pixel = True)

            ## interpolate sun angles to current scene
            #interp_sza = scipy.interpolate.RegularGridInterpolator((x_grid, y_grid), sun_angles[zi, :, :])
            #sza = interp_sza((band_x_mesh, band_y_mesh))
            #interp_saa = scipy.interpolate.RegularGridInterpolator((x_grid, y_grid), sun_angles[ai, :, :])
            #saa = interp_saa((band_x_mesh, band_y_mesh))

            ## interpolate sun angles to current scene - correct y,x order
            interp_sza = scipy.interpolate.RegularGridInterpolator((y_grid, x_grid), sun_angles[zi, :, :],
                                                                  bounds_error = False, fill_value = np.nan)
            sza = interp_sza((band_y_mesh, band_x_mesh))
            interp_saa = scipy.interpolate.RegularGridInterpolator((y_grid, x_grid), sun_angles[ai, :, :],
                                                                  bounds_error = False, fill_value = np.nan)
            saa = interp_saa((band_y_mesh, band_x_mesh))

            # ## warp to target scene
            #if 'region' in dct_sub:
            #    if dct_sub['region'] != dct_prj:
            #    sza = ac.shared.warp_from_source(dct_prj, dct_sub, sza, fill_value = np.nan)
            #    saa = ac.shared.warp_from_source(dct_prj, dct_sub, saa, fill_value = np.nan)

            #if dct_sub != dct_prj:
            if True:
                #sza = ac.shared.warp_from_source(dct_sub, dct_prj, sza, warp_to = warp_to)
                #saa = ac.shared.warp_from_source(dct_sub, dct_prj, saa, warp_to = warp_to)
                #sza = ac.shared.warp_from_source(dct_prj, dct_sub, sza, warp_to = warp_to, source_srs = source_srs)
                #saa = ac.shared.warp_from_source(dct_prj, dct_sub, saa, warp_to = warp_to, source_srs = source_srs)
                #sza = ac.shared.warp_from_source(dct_prj, dct_sub, sza, fill_value = np.nan)
                #saa = ac.shared.warp_from_source(dct_prj, dct_sub, saa, fill_value = np.nan)
                #sza = ac.shared.warp_from_source(dct_sub, dct_prj, sza, warp_to = warp_to, source_srs = source_srs, fill_value = np.nan)
                #saa = ac.shared.warp_from_source(dct_sub, dct_prj, saa, warp_to = warp_to, source_srs = source_srs, fill_value = np.nan)
                sza = ac.shared.warp_from_source(dct, dct_prj, sza, warp_to = warp_to, source_srs = source_srs, fill_value = np.nan)
                saa = ac.shared.warp_from_source(dct, dct_prj, saa, warp_to = warp_to, source_srs = source_srs, fill_value = np.nan)

            ## write data
            gemo.write('sza', sza, replace_nan = True,)
            if verbosity > 1: print('Wrote sza {}'.format(sza.shape))
            sza = None
            gemo.write('saa', saa, replace_nan = True,)
            if verbosity > 1: print('Wrote saa {}'.format(saa.shape))
            ## keep saa as we use it to compute raa
            ##saa = None

        ## extract auxiliary data
        if setu['s2_auxiliary_include']:
            aux_dict = {} ## new dict to store current AUX data
            for aux_source in ['cams', 'ecmwf']:
                aux_variables = [m[0] for m in z['conditions']['meteorology'][aux_source].members()]
                for aux_variable in aux_variables:
                    aux_data = z['conditions']['meteorology'][aux_source][aux_variable][:]
                    aux_key = 'AUX_{}_{}'.format(aux_source.upper(), aux_variable)
                    if aux_key not in gatts: gatts[aux_key] = []
                    #gatts[aux_key] += [v for v in aux_data]
                    gatts[aux_key] += [float(v) for v in aux_data.flatten()]
                    aux_dict[aux_key] = aux_data
                    ## if we want to stack aux data from multiple tiles
                    #if aux_key not in aux_dict:
                    #    aux_dict[aux_key] = aux_data
                    #else:
                    #    ## figure out how to stack aux_data
        ## end extract auxiliary data

        ## delete mask
        mask = None

        ## get band in each resolution group
        bands_10m = [m[0] for m in z['measurements']['reflectance']['r10m'].members() if m[0][0]=='b']
        bands_20m = [m[0] for m in z['measurements']['reflectance']['r20m'].members() if m[0][0]=='b']
        bands_60m = [m[0] for m in z['measurements']['reflectance']['r60m'].members() if m[0][0]=='b']

        ## write TOA bands
        for bi, b in enumerate(rsr_bands):
            #if b != '1': continue ## check b1

            bname = 'B{}'.format(b.zfill(2))
            bname = bname.lower()

            ## get reflectance grid key
            if bname in bands_10m: grid_key = 'r10m'
            elif bname in bands_20m: grid_key = 'r20m'
            elif bname in bands_60m: grid_key = 'r60m'
            else:
                print('{} not found'.format(bname))
                continue
            print('ZARR band name {} and grid_key {}'.format(bname, grid_key))

            ## output geometry
            if setu['output_geometry']:
                ## get detector footprints if B1 or if using per band footprints
                get_footprint = (setu['geometry_fixed_footprint'] & (b == '1')) |\
                                (not setu['geometry_fixed_footprint']) | (setu['geometry_per_band'])
                if get_footprint:
                    ## set up path for current band detector footprint
                    ifile_dfoo = '{}/conditions/mask/detector_footprint/{}/{}'.format(bundle, grid_key, bname)
                    if not url:
                        gdal_file_dfoo = 'ZARR:"{}"'.format(ifile_dfoo)
                    else:
                        gdal_file_dfoo = 'ZARR:"/vsicurl/{}"'.format(ifile_dfoo)

                    ## read data and metadata
                    print('Reading detector footprint for band {}'.format(b))
                    #md, dfoo = ac.shared.read_band(gdal_file_dfoo, sub = sub, warp_to = warp_to, source_srs = source_srs, gdal_meta = True)
                    #md, dfoo = ac.shared.read_band(gdal_file_dfoo, sub = sub, warp_to = None, source_srs = None, gdal_meta = True)

                    ## read full tile dfoo
                    #md, dfoo = ac.shared.read_band(gdal_file_dfoo, sub = None, warp_to = None, source_srs = None, gdal_meta = True)
                    ## read full tile dfoo at 60 metres
                    md, dfoo = ac.shared.read_band(gdal_file_dfoo, sub = None, gdal_meta = True,
                                                   warp_to = warp_to_geom, source_srs = source_srs, )
                    ## find detectors
                    dfoos = np.unique(dfoo)
                    dfoos_names = ['d{}'.format(str(d).zfill(2)) for d in dfoos]
                ## end get detector footprints

                ## compute average geometry if first band
                if (b == '1'):
                    ## compute average view angles grid across the bands (axis 0)
                    viewing_incidence_angles_mean = np.nanmean(viewing_incidence_angles, axis = 0)
                    ## set up mean arrays
                    mean_vza = np.zeros(band_x_mesh.shape) + np.nan
                    mean_vaa = np.zeros(band_x_mesh.shape) + np.nan

                    if plane_fit_geom:
                        ## plane fitting
                        ## run through detectors
                        for di, det_name in enumerate(detectors):
                            if str(det_name) == det_name:
                                det = int(det_name[1:])
                            else:
                                det = int(det_name)
                            valid = np.where(dfoo == det)
                            if len(valid[0]) == 0: continue
                            ## vza
                            normal, d = ac.shared.plane_fit(x_grid_mesh, y_grid_mesh, viewing_incidence_angles_mean[di, zi, :, :], mesh = True)
                            mean_vza[valid] = (-normal[0] * band_x_mesh[valid] - normal[1] * band_y_mesh[valid] + d) / normal[2]
                            ## vaa
                            normal, d = ac.shared.plane_fit(x_grid_mesh, y_grid_mesh, viewing_incidence_angles_mean[di, ai, :, :], mesh = True)
                            mean_vaa[valid] = (-normal[0] * band_x_mesh[valid] - normal[1] * band_y_mesh[valid] + d) / normal[2]
                    else:
                        ## old grid extend method
                        ## run through detectors
                        for di, det_name in enumerate(detectors):
                            if str(det_name) == det_name:
                                det = int(det_name[1:])
                            else:
                                det = int(det_name)
                            det_mask = dfoo == det
                            valid = np.where(det_mask)
                            if len(valid[0]) == 0: continue

                            ## extend the grid
                            vza_ = ac.sentinel2.grid_extend(viewing_incidence_angles_mean[di, zi, :, :], iterations = 1, crop = False)
                            vaa_ = ac.sentinel2.grid_extend(viewing_incidence_angles_mean[di, ai, :, :], iterations = 1, crop = False)

                            ## interpolate
                            mean_vza[valid] = ac.shared.tiles_interp(vza_, xnew+1, ynew+1, smooth = False, fill_nan = True,
                                                              target_mask = det_mask, target_mask_full = False, method='linear')
                            mean_vaa[valid] = ac.shared.tiles_interp(vaa_, xnew+1, ynew+1, smooth = False, fill_nan = True,
                                                              target_mask = det_mask, target_mask_full = False, method='linear')

                    ## warp to target scene
                    #if dct_sub != dct_prj:
                    if True:
                        #mean_vza = ac.shared.warp_from_source(dct_sub, dct_prj, mean_vza, warp_to = warp_to)
                        #mean_vaa = ac.shared.warp_from_source(dct_sub, dct_prj, mean_vaa, warp_to = warp_to)
                        #mean_vza = ac.shared.warp_from_source(dct_prj, dct_sub, mean_vza, warp_to = warp_to, source_srs = source_srs)
                        #mean_vaa = ac.shared.warp_from_source(dct_prj, dct_sub, mean_vaa, warp_to = warp_to, source_srs = source_srs)
                        #mean_vza = ac.shared.warp_from_source(dct_prj, dct_sub, mean_vza, fill_value = np.nan)
                        #mean_vaa = ac.shared.warp_from_source(dct_prj, dct_sub, mean_vaa, fill_value = np.nan)
                        #mean_vza = ac.shared.warp_from_source(dct_sub, dct_prj, mean_vza, warp_to = warp_to, source_srs = source_srs, fill_value = np.nan)
                        #mean_vaa = ac.shared.warp_from_source(dct_sub, dct_prj, mean_vaa, warp_to = warp_to, source_srs = source_srs, fill_value = np.nan)
                        mean_vza = ac.shared.warp_from_source(dct, dct_prj, mean_vza, warp_to = warp_to, source_srs = source_srs, fill_value = np.nan)
                        mean_vaa = ac.shared.warp_from_source(dct, dct_prj, mean_vaa, warp_to = warp_to, source_srs = source_srs, fill_value = np.nan)

                    ## write vza
                    ds = 'vza'
                    gemo.write(ds, mean_vza, replace_nan = True,)
                    if verbosity > 1: print('Wrote {} {}'.format(ds, mean_vza.shape))
                    del mean_vza
                    ## write vaa
                    if setu['s2_write_vaa']:
                        ds = 'vaa'
                        gemo.write(ds, mean_vaa, replace_nan = True,)
                        if verbosity > 1: print('Wrote {} {}'.format(ds, mean_vaa.shape))
                    ## compute raa
                    mean_raa = np.abs(saa - mean_vaa)
                    del mean_vaa
                    tmp = np.where(mean_raa>180)
                    mean_raa[tmp]=np.abs(360 - mean_raa[tmp])
                    ## write raa
                    ds = 'raa'
                    gemo.write(ds, mean_raa, replace_nan = True,)
                    if verbosity > 1: print('Wrote {} {}'.format(ds, mean_raa.shape))
                    del mean_raa, tmp

                    ## write detector footprint
                    if (setu['s2_write_dfoo']):
                        #if dct_sub != dct_prj: ## warp to target scene
                        if True:
                            #dfoo_ = ac.shared.warp_from_source(dct_sub, dct_prj, dfoo, warp_to = warp_to)
                            #dfoo_ = ac.shared.warp_from_source(dct_prj, dct_sub, dfoo, fill_value = np.nan)
                            #dfoo_ = ac.shared.warp_from_source(dct_sub, dct_prj, dfoo, warp_to = warp_to, source_srs = source_srs, fill_value = np.nan)
                            dfoo_ = ac.shared.warp_from_source(dct, dct_prj, dfoo, warp_to = warp_to, source_srs = source_srs, fill_value = np.nan)
                        else:
                            dfoo_ = dfoo * 1
                        ds = 'dfoo'
                        gemo.write(ds, dfoo_, replace_nan = True,)
                        if verbosity > 1: print('Wrote {} {}'.format(ds, dfoo_.shape))
                        del dfoo_
                ## end mean view geometry

                ## start per band geometry if requested
                if setu['geometry_per_band']:
                    print('Computing view geometry for band {}'.format(b))
                    ## empty angle datasets
                    band_vza = np.zeros(band_x_mesh.shape) + np.nan
                    band_vaa = np.zeros(band_x_mesh.shape) + np.nan

                    if plane_fit_geom:
                        ## plane fit
                        ## run through detectors
                        for di, det_name in enumerate(detectors):
                            if str(det_name) == det_name:
                                det = int(det_name[1:])
                            else:
                                det = int(det_name)
                            valid = np.where(dfoo == det)
                            if len(valid[0]) == 0: continue
                            ## vza
                            normal, d = ac.shared.plane_fit(x_grid_mesh, y_grid_mesh, viewing_incidence_angles[bi, di, zi, :, :], mesh = True)
                            band_vza[valid] = (-normal[0] * band_x_mesh[valid] - normal[1] * band_y_mesh[valid] + d) / normal[2]
                            ## vaa
                            normal, d = ac.shared.plane_fit(x_grid_mesh, y_grid_mesh, viewing_incidence_angles[bi, di, ai, :, :], mesh = True)
                            band_vaa[valid] = (-normal[0] * band_x_mesh[valid] - normal[1] * band_y_mesh[valid] + d) / normal[2]
                    else:
                        ## old grid extend method
                        ## run through detectors
                        for di, det_name in enumerate(detectors):
                            if str(det_name) == det_name:
                                det = int(det_name[1:])
                            else:
                                det = int(det_name)
                            det_mask = dfoo == det
                            valid = np.where(det_mask)
                            if len(valid[0]) == 0: continue

                            ## extend the grid
                            vza_ = ac.sentinel2.grid_extend(viewing_incidence_angles[bi, di, zi, :, :], iterations = 1, crop = False)
                            vaa_ = ac.sentinel2.grid_extend(viewing_incidence_angles[bi, di, ai, :, :], iterations = 1, crop = False)
                            ## interpolate
                            band_vza[valid] = ac.shared.tiles_interp(vza_, xnew+1, ynew+1, smooth = False, fill_nan = True,
                                                              target_mask = det_mask, target_mask_full = False, method='linear')
                            band_vaa[valid] = ac.shared.tiles_interp(vaa_, xnew+1, ynew+1, smooth = False, fill_nan = True,
                                                              target_mask = det_mask, target_mask_full = False, method='linear')

                    ## warp to target scene
                    #if dct_sub != dct_prj:
                    if True:
                        #band_vza = ac.shared.warp_from_source(dct_sub, dct_prj, band_vza, warp_to = warp_to)
                        #band_vaa = ac.shared.warp_from_source(dct_sub, dct_prj, band_vaa, warp_to = warp_to)
                        #band_vza = ac.shared.warp_from_source(dct_prj, dct_sub, band_vza, fill_value = np.nan)
                        #band_vaa = ac.shared.warp_from_source(dct_prj, dct_sub, band_vaa, fill_value = np.nan)
                        band_vza = ac.shared.warp_from_source(dct, dct_prj, band_vza, warp_to = warp_to, source_srs = source_srs, fill_value = np.nan)
                        band_vaa = ac.shared.warp_from_source(dct, dct_prj, band_vaa, warp_to = warp_to, source_srs = source_srs, fill_value = np.nan)

                    ## write vza
                    ds = 'vza_{}'.format(rsrd['wave_name'][b])
                    gemo.write(ds, band_vza, replace_nan = True,)
                    if verbosity > 1: print('Wrote {} {}'.format(ds, band_vza.shape))
                    del band_vza

                    ## write vaa
                    if setu['s2_write_vaa']:
                        ds = 'vaa_{}'.format(rsrd['wave_name'][b])
                        gemo.write(ds, band_vaa, replace_nan = True,)
                        if verbosity > 1: print('Wrote {} {}'.format(ds, band_vaa.shape))

                    ## compute raa
                    band_raa = np.abs(saa - band_vaa)
                    del band_vaa
                    tmp = np.where(band_raa>180)
                    band_raa[tmp]=np.abs(360 - band_raa[tmp])

                    ## write raa
                    ds = 'raa_{}'.format(rsrd['wave_name'][b])
                    gemo.write(ds, band_raa, replace_nan = True,)
                    if verbosity > 1: print('Wrote {} {}'.format(ds, band_raa.shape))
                    del band_raa, tmp
                    ## end get view geometry

                    ## write band specific detector footprint
                    if (setu['s2_write_dfoo']) & (setu['s2_write_dfoo_per_band']):
                        if True:
                        #if dct_sub != dct_prj: ## warp to target scene
                        #    #dfoo = ac.shared.warp_from_source(dct_sub, dct_prj, dfoo, warp_to = warp_to)
                        #    dfoo = ac.shared.warp_from_source(dct_prj, dct_sub, dfoo, fill_value = np.nan)
                            dfoo = ac.shared.warp_from_source(dct_sub, dct_prj, dfoo, warp_to = warp_to, source_srs = source_srs, fill_value = np.nan)
                        ds = 'dfoo_{}'.format(rsrd['wave_name'][b])
                        gemo.write(ds, dfoo, replace_nan = True,)
                        if verbosity > 1: print('Wrote {} {}'.format(ds, dfoo.shape))
                    #del dfoo ## don't delete we may be reusing dfoo
            ## end output geometry

            print('Reading rhot for band {}'.format(b))

            ## set up path for current band
            ifile = '{}/measurements/reflectance/{}/{}'.format(bundle, grid_key, bname)
            if not url:
                gdal_file = 'ZARR:"{}"'.format(ifile)
            else:
                gdal_file = 'ZARR:"/vsicurl/{}"'.format(ifile)

            ## read data and metadata
            md, data = ac.shared.read_band(gdal_file, sub = sub, warp_to = warp_to, source_srs = source_srs, gdal_meta = True)

            ## compute mask
            if 'fill_value' in md['_eopf_attrs']:
                data_mask = data == md['_eopf_attrs']['fill_value']
                if setu['s2_dilate_blackfill']:
                    data_mask = scipy.ndimage.binary_dilation(data_mask, iterations = setu['s2_dilate_blackfill_iterations'])

            ## compute reflectance
            ## QV 2025-11-13 these appear to be applied already for current scenes?
            ## QV 2025-11-17 check if dtype is uint16 for scaling
            if data.dtype == 'uint16':
                data = data.astype(np.float32) * md['_eopf_attrs']['scale_factor']
                data += md['_eopf_attrs']['add_offset']

            ## add masks
            data[data_mask] = np.nan
            if (setu['polygon_clip']): data[clip_mask] = np.nan

            ## output dataset and attributes
            ds = 'rhot_{}'.format(rsrd['wave_name'][b])
            ds_att = {'wavelength':rsrd['wave_mu'][b]*1000}

            ## apply gains if provided
            if setu['gains'] & (gains_dict is not None):
                ds_att['toa_gain'] = gains_dict[b]
                data *= ds_att['toa_gain']
                if verbosity > 1: print('Converting bands: Applied TOA gain {} to {}'.format(ds_att['toa_gain'], ds))

            ## apply offsets if provided
            if setu['offsets'] & (offsets_dict is not None):
                ds_att['toa_offset'] = offsets_dict[b]
                data *= ds_att['toa_offset']
                if verbosity > 1: print('Converting bands: Applied TOA offset {} to {}'.format(ds_att['toa_gain'], ds))

            ## write to ms file
            gemo.write(ds, data, replace_nan = True, ds_att = ds_att)
            if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data.shape))
            #stop


        ## update attributes
        gemo.gatts = {k: gatts[k] for k in gatts}
        gemo.gatts_update()
        gemo.close()

        ## reproject AUX data
        ## note these look different to projected GRIB data from SAFE bundle
        ## it is possible those were wrong!
        if (setu['s2_auxiliary_include']) & (setu['s2_auxiliary_project']):
            ## create new aux NetCDF
            if ofile_aux_new:
                ofile_aux = '{}/{}'.format(os.path.dirname(ofile), os.path.basename(ofile).replace('_L1R.nc', '_AUX.nc'))
                gemoa = ac.gem.gem(ofile_aux, new = True)
                gemoa.gatts = {k: gatts[k] for k in gatts}
                gemoa.nc_projection = nc_projection
                ofile_aux_new = False

            ## compute target lon/lat and shape
            lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel=True)
            aux_shape = lon.shape

            ## run through aux products and reproject
            for an in aux_dict:
                if (an.endswith('longitude') | an.endswith('latitude')): continue
                for lkey in ['AUX_CAMS', 'AUX_ECMWF']:
                    if an.startswith(lkey): break

                #lli = np.stack((aux_dict['{}_longitude'.format(lkey)].flatten(),
                #                aux_dict['{}_latitude'.format(lkey)].flatten())).T

                #xi, yi = np.meshgrid(aux_dict['{}_longitude'.format(lkey)],
                #                     aux_dict['{}_latitude'.format(lkey)], indexing='xy')

                ## use bounds_error = False, fill_value = None to extrapolate
                interp_an = scipy.interpolate.RegularGridInterpolator((aux_dict['{}_longitude'.format(lkey)],
                                                                       aux_dict['{}_latitude'.format(lkey)]), aux_dict[an],
                                                                       bounds_error = False, fill_value = None)
                ret = interp_an((lon, lat))
                #ret[mask] = np.nan

                # ## interpolate and fill edges
                # if setu['s2_auxiliary_fill_edge']:
                #     #ret = scipy.interpolate.griddata(lli, aux_dict[an].flatten(), llo)
                #     ret = ac.shared.fillnan(ret.reshape(aux_shape[0], aux_shape[1]))
                #     ret[mask] = np.nan

                ## write
                gemoa.write(an, ret, replace_nan = True)
                if verbosity > 1: print('Wrote {} {}'.format(an, ret.shape))
                del interp_an, ret

        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))

        if setu['limit'] is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
