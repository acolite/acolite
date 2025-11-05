## def l1_convert
## converts Pléiades data to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-02-24
## modifications: 2021-12-31 (QV) new handling of settings
##                2022-01-04 (QV) added netcdf compression
##                2022-11-09 (QV) updates for PNEO processing
##                2023-02-21 (QV) new handling for unprojected and projected data
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-17 (QV) use new gem NetCDF handling
##                2025-01-28 (QV) fix when pan meta is missing
##                2025-01-30 (QV) moved polygon limit
##                2025-02-02 (QV) removed percentiles
##                2025-02-04 (QV) improved settings handling
##                2025-02-07 (QV) added pleiades_force_metadata_geolocation
##                2025-02-10 (QV) cleaned up settings use
##                2025-09-15 (QV) fixed processing of Airbus "REFLECTANCE" data
##                                updated data reading and added half pixel to nc_projection

def l1_convert(inputfile, output = None, settings = None):

    import os, copy
    import dateutil.parser, time
    import numpy as np
    import acolite as ac
    import scipy.ndimage

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

    ofiles = []
    for bundle in inputfile:
        ofile = None

        ## check if we are dealing with Pléiades image bundle
        try:
            ifiles,mfiles,pifiles,pmfiles = ac.pleiades.bundle_test(bundle, listpan=True)
        except:
            print('Bundle {} not recognised.'.format(bundle))
            continue

        ## read meta data
        if verbosity > 1: print('Importing metadata from {}'.format(bundle))
        mfiles_set = set(mfiles)
        if len(mfiles_set)>1:
            print('Multiple metadata files found')
            return()

        sub = None
        pansub = None
        new = True
        new_pan = True
        ofile = None

        t0 = time.time()

        ## read metadata only once
        pmeta = None
        for mfile in mfiles_set: meta = ac.pleiades.metadata_parse(mfile)
        for pmfile in set(pmfiles):
            if pmfile == '': continue
            pmeta = ac.pleiades.metadata_parse(pmfile, pan=True)
        sensor = meta['sensor']

        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults

        verbosity = setu['verbosity']
        if output is None: output = setu['output']
        if output is None: output = os.path.dirname(bundle)

        ## set resolution and band names
        if meta['sensor'] in ['PHR1A', 'PHR1B']:
            btags = {'Blue':'B0', 'Green':'B1', 'Red':'B2', 'NIR':'B3', 'Pan':'P'}
            ms_resolution = 2
            p_resolution = 0.5
        elif meta['sensor'] in ['PNEO3', 'PNEO4']:
            btags = {'BlueCoastal':'DB', 'Blue':'B', 'Green':'G', 'Red':'R', 'RedEdge':'RE', 'NIR':'NIR', 'PAN':'PAN'}
            ms_resolution = 1.32
            p_resolution = 0.33

        ## parse metadata
        dtime = dateutil.parser.parse(meta['isotime'])
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()

        ## read rsr
        sensor = meta['sensor']
        rsrd = ac.shared.rsr_dict(sensor)[sensor]

        ## get F0 - not stricty necessary if using provided F0
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        f0_b = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data']*10, rsrd['rsr'])

        gatts = {'sensor':meta['sensor'], 'satellite':meta['satellite'],
                 'satellite_sensor':'{}_{}'.format(meta['satellite'], meta['sensor']),
                 'isodate':isodate,
                 'sza':meta['sza'], 'vza':meta['vza'], 'saa':meta['saa'],'vaa':meta['vaa'],
                 'raa':meta['raa'], 'se_distance': se_distance,
                 'mus': np.cos(meta['sza']*(np.pi/180.)),
                 'acolite_file_type': 'L1R'}

        ## add band info to gatts
        for b in rsrd['rsr_bands']:
            gatts['{}_wave'.format(b)] = rsrd['wave_nm'][b]
            gatts['{}_name'.format(b)] = rsrd['wave_name'][b]
            gatts['{}_f0'.format(b)] = f0_b[b]

        stime = dateutil.parser.parse(gatts['isodate'])

        ## set up oname (without directory or file type) and ofile (with directory and file type)
        oname = '{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        oname = '{}_{}'.format(oname, gatts['acolite_file_type'])
        ofile = '{}/{}.nc'.format(output, oname)
        pofile = '{}/{}_pan.nc'.format(output, oname)
        gatts['oname'] = oname
        gatts['ofile'] = ofile
        gatts['ofile_pan'] = pofile

        ## test scene
        if setu['limit'] is not None:
            out_scene = ac.pleiades.geo.test_coverage(meta, setu['limit'], verbose=verbosity>2)
            if out_scene:
                print('Provided limit {} not covered by scene'.format(setu['limit']))
                continue

        dct = None
        ## find projection info based on first file
        ## of pleiades_force_metadata_geolocation is set the "old method" will be used
        if not setu['pleiades_force_metadata_geolocation']:
            ifile = ifiles[0]
            try:
                print('Reading projection from {}'.format(ifile))
                dct = ac.shared.projection_read(ifile, use_world_transform = False)
            except:
                print('Could not determine projection from {}'.format(ifile))

            ## if we can read the projection, read other tiles and update dct
            if dct is not None:
                for ifile in ifiles[1:]:
                    ## read projection
                    print('Reading projection from {}'.format(ifile))
                    dct_tile = ac.shared.projection_read(ifile, use_world_transform = False)
                    ## compute new ranges
                    dct['xrange'] = min(dct['xrange'][0], dct_tile['xrange'][0]),\
                                        max(dct['xrange'][1], dct_tile['xrange'][1])
                    dct['yrange'] = max(dct['yrange'][0], dct_tile['yrange'][0]),\
                                        min(dct['yrange'][1], dct_tile['yrange'][1])
                ## update dimensions
                dct['xdim'] = np.round((dct['xrange'][1] - dct['xrange'][0]) / dct['pixel_size'][0]).astype(int)
                dct['ydim'] = np.round((dct['yrange'][1] - dct['yrange'][0]) / dct['pixel_size'][1]).astype(int)
                dct['dimensions'] = (dct['xdim'], dct['ydim'])
                print(dct)

        rpc_use = False ## if projected and RPC are given this will lead to offsets
        new_method, reproject = False, False
        ## set up reprojection
        if (setu['reproject_inputfile_force']) | (setu['reproject_inputfile'] & (dct is None)):
            ## get image extent from metadata
            vlons = [meta['VERTICES'][v]['LON'] for v in meta['VERTICES']]
            vlats = [meta['VERTICES'][v]['LAT'] for v in meta['VERTICES']]

            ## set up ms projection
            if setu['limit'] is None:
                vlimit = [np.nanmin(vlats), np.nanmin(vlons), np.nanmax(vlats), np.nanmax(vlons)]
                dct, nc_projection, warp_to = ac.shared.projection_setup(vlimit, ms_resolution, res_method = setu['warp_resampling_method'])
            else:
                dct, nc_projection, warp_to = ac.shared.projection_setup(setu['limit'], ms_resolution, res_method = setu['warp_resampling_method'])
            reproject = True
            ## update gatts
#             gatts['scene_xrange'] = dct['xrange']
#             gatts['scene_yrange'] = dct['yrange']
#             gatts['scene_proj4_string'] = dct['proj4_string']
#             gatts['scene_pixel_size'] = dct['pixel_size']
#             gatts['scene_dims'] = dct['dimensions']
#             if 'zone' in dct: gatts['scene_zone'] = dct['zone']

        add_half_pixel = True
        targetAlignedPixels = False
        if (dct is not None):
            new_method = True
            res_method = 'average'
            ## check limit
            if (setu['limit'] is not None) & (reproject is False):
                dct_sub = ac.shared.projection_sub(dct, setu['limit'], four_corners = True)
                if dct_sub['out_lon']:
                    if verbosity > 1: print('Longitude limits outside {}'.format(bundle))
                    continue
                if dct_sub['out_lat']:
                    if verbosity > 1: print('Latitude limits outside {}'.format(bundle))
                    continue
                dct = {k: dct_sub[k] for k in dct_sub}

            ## MS data
            nc_projection = ac.shared.projection_netcdf(dct, add_half_pixel = add_half_pixel)
            xyr = [min(dct['xrange']),min(dct['yrange']),
                   max(dct['xrange']),max(dct['yrange']), dct['projection']]
            warp_to = (dct['projection'], xyr, dct['pixel_size'][0],dct['pixel_size'][1], res_method)
            #print(warp_to)
            #print(nc_projection)

            ## PAN data
            dct_pan = {k:dct[k] for k in dct}
            dct_pan['pixel_size'] = (dct['pixel_size'][0] * (p_resolution/ms_resolution), dct['pixel_size'][1] * (p_resolution/ms_resolution))
            dct_pan['dimensions'] = np.round((dct_pan['yrange'][1]-dct_pan['yrange'][0]) / dct_pan['pixel_size'][1],0).astype(int), \
                                    np.round((dct_pan['xrange'][1]-dct_pan['xrange'][0]) / dct_pan['pixel_size'][0],0).astype(int)
            dct_pan['ydim'] = dct_pan['dimensions'][0]
            dct_pan['xdim'] = dct_pan['dimensions'][1]
            nc_projection_pan = ac.shared.projection_netcdf(dct_pan, add_half_pixel = add_half_pixel)
            xyr_pan = [min(dct_pan['xrange']), min(dct_pan['yrange']),
                       max(dct_pan['xrange']), max(dct_pan['yrange']), dct_pan['projection']]
            warp_to_pan = (dct_pan['projection'], xyr_pan, dct_pan['pixel_size'][0],dct_pan['pixel_size'][1], res_method)

            ## set up RPC DEM
            rpc_dem = None
            if setu['reproject_inputfile_dem']: rpc_dem = ac.dem.copernicus_dem_rpc(dct, output=output)

        ## convert using "old" method
        ## run through tiles and subset output data
        if (dct is None) & (new_method is False):
            print('Pléiades old method')
            t = time.process_time()
            gatts['scene_dims'] = int(meta['NROWS']), int(meta['NCOLS'])

            ## check if we need to find crop position
            if setu['limit'] is not None:
                if len(setu['limit']) == 4:
                    ncols = int(meta['NCOLS'])
                    nrows = int(meta['NROWS'])
                    sub = ac.pleiades.geo.crop(meta, setu['limit'])
                    if pmeta is not None:
                        pansub = ac.pleiades.geo.crop(pmeta, setu['limit'])
                        ## QV 2022-11-09
                        ## force pan sub dimensions to match ms data
                        ## to be improved!
                        pansub[2] = sub[2]*4
                        pansub[3] = sub[3]*4

            if sub is not None:
                sub = [int(s) for s in sub]
                cropname = '_'.join([str(i) for i in sub])

            if sub is None:
                dims = int(meta['NROWS']), int(meta['NCOLS'])
                gatts['global_dims'] = dims
            else:
                dims = sub[3], sub[2]
                gatts['global_dims'] = sub[3], sub[2]

            new = True
            if new:
                gemo = ac.gem.gem(ofile, new = True)
                gemo.gatts = {k: gatts[k] for k in gatts}
                datasets = gemo.datasets

            ## write lat/lon
            if (setu['output_geolocation']):
                if ('lat' not in datasets) or ('lon' not in datasets):
                    if verbosity > 1: print('Writing geolocation lon/lat')
                    lon, lat = ac.pleiades.geo.ll(meta, sub=sub)
                    gemo.write('lon', lon)
                    #ac.output.nc_write(ofile, 'lon', lon, attributes=gatts, new=new)
                    lon = None
                    if verbosity > 1: print('Wrote lon')
                    #ac.output.nc_write(ofile, 'lat', lat)
                    gemo.write('lat', lat)
                    lat = None
                    if verbosity > 1: print('Wrote lat')
                    new = False

            ## run through image tiles
            tile_id = 0
            tile_dims = meta['tiles_nrows'], meta['tiles_ncols']

            for r_tile in range(meta['ntiles_R']):
                for c_tile in range(meta['ntiles_C']):
                    tile_id+=1
                    tile_name = 'R{}C{}'.format(r_tile+1,c_tile+1)
                    print(tile_id, tile_name)

                    ## tile offsets
                    tile_row_off = r_tile * tile_dims[0] # Y
                    tile_col_off = c_tile * tile_dims[1] # X

                    ## subsetting tiled image
                    ## not tested for ROI across tiles!
                    sub_tile = None
                    pansub_tile = None
                    if sub is not None:
                        sub_tile = [sub[0]-tile_col_off, sub[1]-tile_row_off,sub[2], sub[3]]

                        if sub_tile[0] < 0: continue
                        if sub_tile[1] < 0: continue
                        if sub_tile[0] > tile_dims[0]: continue
                        if sub_tile[1] > tile_dims[1]: continue
                        pansub_tile = [pansub[0]-(tile_col_off*4),
                                       pansub[1]-(tile_row_off*4),
                                       pansub[2], pansub[3]]

                    ifile = None
                    for it,tfile in enumerate(ifiles):
                        if tile_name not in tfile: continue
                        tbn = os.path.basename(tfile)
                        mfile=mfiles[it]
                        pifile=pifiles[it]
                        pmfile=pmfiles[it]
                        ## track PNEO RGB and NED tiles
                        if ('PNEO' in meta['sensor']):
                            if '_NED_' in tbn: ifile_ned=ifiles[it]
                            if '_RGB_' in tbn: ifile=ifiles[it]
                        else:
                            ifile=ifiles[it]

                    dct = None
                    nc_projection = None
                    update_projection = True

                    dct_pan = None
                    nc_projection_pan = None
                    update_projection_pan = True

                    ## read in TOA reflectances
                    for b in rsrd['rsr_bands']:
                        pan = False
                        if btags[b] in meta['BAND_INFO']:
                            bd = {k:meta['BAND_INFO'][btags[b]][k] for k in meta['BAND_INFO'][btags[b]]}
                            #idx = 1+bd['band_index']
                            idx = 0 + bd['band_index']
                            ## read data
                            ifile_ = '{}'.format(ifile)

                            ## Red Green Blue PNEO bands in RGB file
                            ## NIR, RedEdge, CoastalBlue in NED file
                            if ('PNEO' in meta['sensor']) & (b in ['BlueCoastal', 'RedEdge', 'NIR']):
                                ifile_ = '{}'.format(ifile_ned)

                            print('Reading band {} from {}'.format(b, ifile_))
                            data = ac.shared.read_band(ifile_, idx=idx, sub=sub_tile)
                            if update_projection:
                                try:
                                    dct = ac.shared.projection_read(ifile_, use_world_transform = False)
                                    #nc_projection = ac.shared.projection_netcdf(dct, add_half_pixel = add_half_pixel)
                                    if sub is not None:
                                        dct = ac.shared.projection_sub_dct(dct, sub)
                                    nc_projection = ac.shared.projection_netcdf(dct, add_half_pixel = add_half_pixel)
                                except BaseException as err:
                                    print("Could not determine projection from {} error {}".format(ifile_, type(err)))
                                    pass
                        else:
                            if pmeta is None: continue
                            if setu['pleiades_skip_pan']: continue
                            pan = True
                            if sub is None:
                                pandims = int(pmeta['NROWS']), int(pmeta['NCOLS'])
                            else:
                                pandims = pansub[3], pansub[2]

                            print('Reading band {} from {}'.format(b, pifile))
                            ## read data
                            data = ac.shared.read_band(pifile, idx=1, sub=pansub_tile)
                            try:
                                dct_pan = ac.shared.projection_read(pifile, use_world_transform = False)
                                #c_projection_pan = ac.shared.projection_netcdf(dct_pan, add_half_pixel = add_half_pixel)
                                if pansub is not None:
                                    dct_pan = ac.shared.projection_sub_dct(dct_pan, pansub)
                                nc_projection_pan = ac.shared.projection_netcdf(dct_pan, add_half_pixel = add_half_pixel)
                            except BaseException as err:
                                print("Could not determine projection from {} error {}".format(pifile, type(err)))
                                pass

                        nodata = data == np.uint16(meta['NODATA'])
                        nodata2 = data == 1
                        saturated = data >= np.uint16(meta['SATURATED']) -1

                        print(idx, b, btags[b])
                        data = data.astype(np.float32)
                        if (meta['RADIOMETRIC_PROCESSING'] == 'RADIANCE') | (meta['RADIOMETRIC_PROCESSING'] == 'BASIC'):
                            data *= (1./bd['radiance_gain'])
                            data += (bd['radiance_bias'])
                            data *= (np.pi * gatts['se_distance']**2) / (bd['F0'] * gatts['mus'])
                        elif (meta['RADIOMETRIC_PROCESSING'] == 'LINEAR_STRETCH'):
                            print('Warning linear stretch data')
                            data *= (1./bd['radiance_gain'])
                            data += (bd['radiance_bias'])
                            data *= (np.pi * gatts['se_distance']**2) / (bd['F0'] * gatts['mus'])
                        elif (meta['RADIOMETRIC_PROCESSING'] == 'REFLECTANCE'):
                            ## convert to reflectance (with AIRBUS A/C)
                            data /= bd['reflectance_gain']
                            data += bd['reflectance_bias']
                            ## convert to LTOA
                            data /= bd['radiance_gain']
                            data += (bd['radiance_bias'])
                            ## and convert to RHOT
                            data *= (np.pi * gatts['se_distance']**2) / (bd['F0'] * gatts['mus'])
                        else:
                            print("{} RADIOMETRIC_PROCESSING not recognised".format(metadata['RADIOMETRIC_PROCESSING']))
                            continue

                        data[nodata] = np.nan
                        data[nodata2] = np.nan
                        data[saturated] = np.nan

                        ds = 'rhot_{}'.format(rsrd['wave_name'][b])
                        ds_att = {'wavelength':rsrd['wave_nm'][b]}

                        ## QV 2022-11-09
                        ## nc_projection does not match when using crop
                        #if sub is not None:
                        #    nc_projection = None
                        #    nc_projection_pan = None

                        if pan:
                            cur_shape = data.shape
                            if cur_shape != pandims:
                                data_full = np.zeros(pandims)+np.nan
                                data_full[tile_row_off*4:tile_row_off*4+cur_shape[0],
                                          tile_col_off*4:tile_col_off*4+cur_shape[1]] = data
                            else:
                                data_full = data * 1.0

                            ## write to netcdf file
                            if new_pan:
                                gemop = ac.gem.gem(pofile, new = True)
                                gemop.gatts = {k: gatts[k] for k in gatts}
                                gemop.nc_projection = nc_projection_pan

                            if (update_projection_pan) & (nc_projection_pan is not None):
                                print('Updating PAN projection')
                                gemop.nc_projection = nc_projection_pan
                            gemop.write(ds, data_full, ds_att = ds_att, replace_nan=True, update_projection = update_projection_pan)

                            #ac.output.nc_write(pofile, ds, data_full, replace_nan=True, #attributes=gatts,
                            #                    new=new_pan, dataset_attributes = ds_att,
                            #                    nc_projection=nc_projection_pan, update_projection=update_projection_pan)
                            data_full = None
                            update_projection_pan = False
                            new_pan = False

                            ## mask data before zooming
                            dmin = np.nanmin(data)
                            data[np.isnan(data)] = 0
                            data = scipy.ndimage.zoom(data, 0.25, order=1)
                            data[data<dmin] = np.nan
                            data[data==dmin] = np.nan

                        cur_shape = data.shape
                        print(dims)
                        print(cur_shape)
                        print(tile_row_off,tile_col_off)
                        if cur_shape != dims:
                            data_full = np.zeros(dims)+np.nan
                            data_full[tile_row_off:tile_row_off+cur_shape[0],
                                      tile_col_off:tile_col_off+cur_shape[1]] = data
                        else:
                            data_full = data * 1.0
                        data = None

                        ## update lat/lon if we were able to compute nc_projection
                        if (update_projection) & (nc_projection is not None):
                            if verbosity > 1: print('Writing geolocation lon/lat')
                            lon, lat = ac.shared.projection_geo(dct, add_half_pixel = add_half_pixel)
                            print('Updating MS projection')
                            gemo.nc_projection = nc_projection
                            print(lon.shape)
                            gemo.write('lon', lon, update_projection = True)
                            lon = None
                            update_projection = False

                            if verbosity > 1: print('Wrote lon')
                            print(lat.shape)
                            gemo.write('lat', lat)
                            lat = None
                            if verbosity > 1: print('Wrote lat')

                        ## write to netcdf file
                        gemo.write(ds, data_full, ds_att = ds_att, replace_nan=True)
                        update_projection = False
                        new = False
                        if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data_full.shape))
            ### end old method
        else:
            print('Pléiades new method')
            if new:
                gemo = ac.gem.gem(ofile, new = True)
                gemo.gatts = {k: gatts[k] for k in gatts}
                gemo.nc_projection = nc_projection
                datasets = gemo.datasets

            ## run through image tiles
            t = time.process_time()

            ## write lat/lon
            if (setu['output_geolocation']):
                if verbosity > 1: print('Writing geolocation lon/lat')
                lon, lat = ac.shared.projection_geo(dct, add_half_pixel = add_half_pixel)
                gemo.write('lon', lon)
                if verbosity > 1: print('Wrote lon ({})'.format(lon.shape))
                lon = None
                gemo.write('lat', lat)
                if verbosity > 1: print('Wrote lat ({})'.format(lat.shape))
                lat = None
                #new=False

            ## read in TOA reflectances
            for b in rsrd['rsr_bands']:
                print(b)
                bd = None
                if b in ['Pan', 'PAN']:
                    pan = True
                    if pmeta is None: continue
                    if setu['pleiades_skip_pan']: continue
                    tfiles = [f for f in pifiles]
                    bd = {k:pmeta['BAND_INFO'][btags[b]][k] for k in pmeta['BAND_INFO'][btags[b]]}
                    idx = 1
                else:
                    pan = False
                    if ('PNEO' in meta['sensor']):
                        tfiles = [f for f in ifiles if '_RGB_' in os.path.basename(f)]
                        tfiles.sort()
                        tfiles_ned = [f for f in ifiles if '_NED_' in os.path.basename(f)]
                        tfiles_ned.sort()
                    else:
                        tfiles = [f for f in ifiles]
                    bd = {k:meta['BAND_INFO'][btags[b]][k] for k in meta['BAND_INFO'][btags[b]]}
                    idx = 0 + bd['band_index']
                if bd is None:
                    print('{} not found in BAND_INFO'.format(b))
                    continue

                its = 1
                if (pan) & (setu['pleiades_skip_pan'] is False): its = 2

                for ii in range(its):
                    if (ii == 1):
                        if pan is False: continue

                    ## run through tiles
                    for it,tfile in enumerate(tfiles):
                        tfile_ = '{}'.format(tfile)
                        ## Red Green Blue PNEO bands in RGB file
                        ## NIR, RedEdge, CoastalBlue in NED file
                        if ('PNEO' in meta['sensor']) & (b in ['BlueCoastal', 'RedEdge', 'NIR']):
                            tfile_ = '{}'.format(tfiles_ned[it])
                        print('Reading band {} from {}'.format(b, tfile_))

                        if (ii == 0):
                            data_in = ac.shared.read_band(tfile_, idx = idx, sub = sub, warp_to = warp_to,
                                rpc_dem = rpc_dem, rpc_use = rpc_use, targetAlignedPixels = targetAlignedPixels)
                        else:
                            data_in = ac.shared.read_band(tfile_, idx = idx, sub = sub, warp_to = warp_to_pan,
                                rpc_dem = rpc_dem, rpc_use = rpc_use, targetAlignedPixels = targetAlignedPixels)
                        nodata = data_in == np.uint16(meta['NODATA'])
                        nodata2 = data_in == 1
                        saturated = data_in >= np.uint16(meta['SATURATED']) -1

                        data_in = data_in.astype(np.float32)
                        data_in[nodata] = np.nan
                        data_in[nodata2] = np.nan
                        data_in[saturated] = np.nan

                        if it == 0:
                            data = data_in
                        else:
                            #dsub = np.where((data_in>=1) & (data < 1))
                            dsub = np.where((data_in>0) & (data == 0))
                            data[dsub] = data_in[dsub]
                    ## end run through tiles

                    ## convert to reflectance
                    if (meta['RADIOMETRIC_PROCESSING'] == 'RADIANCE') | (meta['RADIOMETRIC_PROCESSING'] == 'BASIC'):
                        data *= (1./bd['radiance_gain'])
                        data += (bd['radiance_bias'])
                        data *= (np.pi * gatts['se_distance']**2) / (bd['F0'] * gatts['mus'])
                    elif (meta['RADIOMETRIC_PROCESSING'] == 'LINEAR_STRETCH'):
                        print('Warning linear stretch data')
                        data *= (1./bd['radiance_gain'])
                        data += (bd['radiance_bias'])
                        data *= (np.pi * gatts['se_distance']**2) / (bd['F0'] * gatts['mus'])
                    elif (meta['RADIOMETRIC_PROCESSING'] == 'REFLECTANCE'):
                        print("REFLECTANCE RADIOMETRIC_PROCESSING, converting to LTOA and RHOT")
                        ## convert to reflectance (with AIRBUS A/C)
                        data /= bd['reflectance_gain']
                        data += bd['reflectance_bias']
                        ## convert to LTOA
                        data /= bd['radiance_gain']
                        data += (bd['radiance_bias'])
                        ## and convert to RHOT
                        data *= (np.pi * gatts['se_distance']**2) / (bd['F0'] * gatts['mus'])
                    else:
                        print("{} RADIOMETRIC_PROCESSING not recognised".format(metadata['RADIOMETRIC_PROCESSING']))
                        continue

                    ds = 'rhot_{}'.format(rsrd['wave_name'][b])
                    ds_att = {'wavelength':rsrd['wave_nm'][b]}

                    ## write to netcdf file
                    if (ii == 0):
                        gemo.write(ds, data, ds_att = ds_att)
                    else: ## second iteration is pan at native resolution
                        ## write to netcdf file
                        gemop = ac.gem.gem(pofile, new = True)
                        gemop.gatts = {k: gatts[k] for k in gatts}
                        gemop.nc_projection = nc_projection_pan
                        gemop.write(ds, data, ds_att = ds_att)

                    if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data.shape))
        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))
        gemo.close()
        try:
            gemop.close()
        except:
            pass
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
