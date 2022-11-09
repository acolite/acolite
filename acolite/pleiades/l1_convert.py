## def l1_convert
## converts Pléiades data to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-02-24
## modifications: 2021-12-31 (QV) new handling of settings
##                2022-01-04 (QV) added netcdf compression
##                2022-11-09 (QV) updates for PNEO processing

def l1_convert(inputfile, output = None, settings = {},
                limit = None, sub = None,
                poly = None,
                output_geolocation = True,
                skip_pan = False,
                percentiles_compute = True,
                percentiles = (0,1,5,10,25,50,75,90,95,99,100),
                verbosity = 5, vname = ''):

    import os
    import dateutil.parser, time
    import numpy as np
    import acolite as ac
    import scipy.ndimage

    if 'verbosity' in settings: verbosity = settings['verbosity']

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

        ## merge sensor specific settings
        setu = ac.acolite.settings.parse(meta['sensor'], settings=settings)
        verbosity = setu['verbosity']
        if output is None: output = setu['output']
        limit = setu['limit']
        poly = setu['polygon']
        output_geolocation = setu['output_geolocation']
        skip_pan = setu['pleiades_skip_pan']

        ## check if ROI polygon is given
        clip, clip_mask = False, None
        if poly is not None:
            if os.path.exists(poly):
                try:
                    limit = ac.shared.polygon_limit(poly)
                    print('Using limit from polygon envelope: {}'.format(limit))
                    clip = True
                except:
                    print('Failed to import polygon {}'.format(poly))

        if limit is not None:
            out_scene = ac.pleiades.geo.test_coverage(meta, limit, verbose=verbosity>2)
            if out_scene:
                print('Provided limit {} not covered by scene'.format(limit))
                continue

        ## align metadata tag names with RSR names
        if 'PNEO' in meta['sensor']:
            btags = {'BlueCoastal':'DB', 'Blue':'B', 'Green':'G', 'Red':'R', 'RedEdge':'RE', 'NIR':'NIR', 'PAN':'PAN'}
        else:
            btags = {'Blue':'B0', 'Green':'B1', 'Red':'B2', 'NIR':'B3', 'Pan':'P'}

        dtime = dateutil.parser.parse(meta['isotime'])
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()

        ## read rsr
        rsrf = ac.path+'/data/RSR/{}.txt'.format(meta['sensor'])
        rsr, rsr_bands = ac.shared.rsr_read(rsrf)
        waves = np.arange(250, 2500)/1000
        waves_mu = ac.shared.rsr_convolute_dict(waves, waves, rsr)
        waves_names = {'{}'.format(b):'{:.0f}'.format(waves_mu[b]*1000) for b in waves_mu}

        ## get F0 - not stricty necessary if using provided F0
        f0 = ac.shared.f0_get()
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data'])*10, rsr)

        gatts = {'sensor':meta['sensor'], 'satellite':meta['satellite'],
                 'satellite_sensor':'{}_{}'.format(meta['satellite'], meta['sensor']),
                 'isodate':isodate,
                 'sza':meta['sza'], 'vza':meta['vza'], 'raa':meta['raa'], 'se_distance': se_distance,
                 'mus': np.cos(meta['sza']*(np.pi/180.)),
                 'acolite_file_type': 'L1R'}

        stime = dateutil.parser.parse(gatts['isodate'])
        oname = '{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'))
        if vname != '': oname+='_{}'.format(vname)

        ## output file information
        ofile = '{}/{}_L1R.nc'.format(output, oname)
        pofile = '{}/{}_L1R_pan.nc'.format(output, oname)
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## add band info to gatts
        for b in rsr_bands:
            gatts['{}_wave'.format(b)] = waves_mu[b]*1000
            gatts['{}_name'.format(b)] = waves_names[b]
            gatts['{}_f0'.format(b)] = f0_b[b]

        gatts['scene_dims'] = int(meta['NROWS']), int(meta['NCOLS'])

        ## check if we need to find crop position
        if limit is not None:
            if len(limit) == 4:
                ncols = int(meta['NCOLS'])
                nrows = int(meta['NROWS'])
                sub = ac.pleiades.geo.crop(meta, limit)
                if pmeta is not None:
                    pansub = ac.pleiades.geo.crop(pmeta, limit)
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
        ## write lat/lon
        if (output_geolocation):
            if (os.path.exists(ofile) & (not new)):
                datasets = ac.shared.nc_datasets(ofile)
            else:
                datasets = []
            if ('lat' not in datasets) or ('lon' not in datasets):
                if verbosity > 1: print('Writing geolocation lon/lat')
                lon, lat = ac.pleiades.geo.ll(meta, sub=sub)
                ac.output.nc_write(ofile, 'lon', lon, attributes=gatts, new=new, double=True,
                                    netcdf_compression=setu['netcdf_compression'],
                                    netcdf_compression_level=setu['netcdf_compression_level'])
                lon = None
                if verbosity > 1: print('Wrote lon')
                ac.output.nc_write(ofile, 'lat', lat, double=True,
                                    netcdf_compression=setu['netcdf_compression'],
                                    netcdf_compression_level=setu['netcdf_compression_level'])
                lat = None
                if verbosity > 1: print('Wrote lat')
                new = False

        ## run through image tiles
        t = time.process_time()
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
                for b in rsr_bands:
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
                                dct = ac.shared.projection_read(ifile_)
                                nc_projection = ac.shared.projection_netcdf(dct, add_half_pixel=False)
                                if sub is not None:
                                    dct = ac.shared.projection_sub_dct(dct, sub)
                                    nc_projection = ac.shared.projection_netcdf(dct, add_half_pixel=False)
                            except BaseException as err:
                                print(f"Could not determine projection from {ifile_=} error {err=}, {ifile_=}, {type(err)=}")
                                pass
                    else:
                        if pmeta is None: continue
                        if skip_pan: continue
                        pan = True
                        if sub is None:
                            pandims = int(pmeta['NROWS']), int(pmeta['NCOLS'])
                        else:
                            pandims = pansub[3], pansub[2]

                        print('Reading band {} from {}'.format(b, pifile))
                        ## read data
                        data = ac.shared.read_band(pifile, idx=1, sub=pansub_tile)
                        try:
                            dct_pan = ac.shared.projection_read(pifile)
                            nc_projection_pan = ac.shared.projection_netcdf(dct_pan, add_half_pixel=False)
                            if pansub is not None:
                                dct_pan = ac.shared.projection_sub_dct(dct_pan, pansub)
                                nc_projection_pan = ac.shared.projection_netcdf(dct_pan, add_half_pixel=False)
                        except BaseException as err:
                            print(f"Could not determine projection from {pifile=} error {err=}, {pifile=}, {type(err)=}")
                            pass

                    nodata = data == np.uint16(meta['NODATA'])
                    nodata2 = data == 1
                    print(idx, b, btags[b])
                    data = data.astype(np.float32)
                    if (meta['RADIOMETRIC_PROCESSING'] == 'RADIANCE') | (meta['RADIOMETRIC_PROCESSING'] == 'BASIC'):
                        #data *= (1./meta['BAND_INFO'][btags[b]]['radiance_gain'])
                        #data += (meta['BAND_INFO'][btags[b]]['radiance_bias'])
                        #data *= (np.pi * gatts['se_distance']**2) / (meta['BAND_INFO'][btags[b]]['F0'] * gatts['mus'])
                        data *= (1./bd['radiance_gain'])
                        data += (bd['radiance_bias'])
                        data *= (np.pi * gatts['se_distance']**2) / (bd['F0'] * gatts['mus'])
                    elif (meta['RADIOMETRIC_PROCESSING'] == 'LINEAR_STRETCH'):
                        print('Warning linear stretch data')
                        #data *= (1./meta['BAND_INFO'][btags[b]]['radiance_gain'])
                        #data += (meta['BAND_INFO'][btags[b]]['radiance_bias'])
                        #data *= (np.pi * gatts['se_distance']**2) / (meta['BAND_INFO'][btags[b]]['F0'] * gatts['mus'])
                        data *= (1./bd['radiance_gain'])
                        data += (bd['radiance_bias'])
                        data *= (np.pi * gatts['se_distance']**2) / (bd['F0'] * gatts['mus'])
                    elif (meta['RADIOMETRIC_PROCESSING'] == 'REFLECTANCE'):
                        #data /= meta['BAND_INFO'][btags[b]]['reflectance_gain']
                        #data += meta['BAND_INFO'][btags[b]]['reflectance_bias']
                        data /= bd['reflectance_gain']
                        data += bd['reflectance_bias']
                        data /= gatts['mus']
                    else:
                        print("{} RADIOMETRIC_PROCESSING not recognised".format(metadata['RADIOMETRIC_PROCESSING']))
                        continue

                    data[nodata] = np.nan
                    data[nodata2] = np.nan

                    ds = 'rhot_{}'.format(waves_names[b])
                    ds_att = {'wavelength':waves_mu[b]*1000}
                    if percentiles_compute:
                        ds_att['percentiles'] = percentiles
                        ds_att['percentiles_data'] = np.nanpercentile(data, percentiles)

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
                        ac.output.nc_write(pofile, ds, data_full, replace_nan=True, #attributes=gatts,
                                            new=new_pan, dataset_attributes = ds_att,
                                            nc_projection=nc_projection_pan, update_projection=update_projection_pan,
                                            netcdf_compression=setu['netcdf_compression'],
                                            netcdf_compression_level=setu['netcdf_compression_level'],
                                            netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
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
                        lon, lat = ac.shared.projection_geo(dct, add_half_pixel=False)
                        print(lon.shape)
                        ac.output.nc_write(ofile, 'lon', lon, double=True,
                                            nc_projection=nc_projection, update_projection=update_projection,
                                            netcdf_compression=setu['netcdf_compression'],
                                            netcdf_compression_level=setu['netcdf_compression_level'])
                        lon = None
                        update_projection = False

                        if verbosity > 1: print('Wrote lon')
                        print(lat.shape)
                        ac.output.nc_write(ofile, 'lat', lat, double=True,
                                            nc_projection=nc_projection, update_projection=update_projection,
                                            netcdf_compression=setu['netcdf_compression'],
                                            netcdf_compression_level=setu['netcdf_compression_level'])
                        lat = None
                        if verbosity > 1: print('Wrote lat')
                        new=False

                    ## write to netcdf file
                    ac.output.nc_write(ofile, ds, data_full, replace_nan=True, attributes=gatts, new=new, dataset_attributes = ds_att,
                                        nc_projection=nc_projection, update_projection=update_projection,
                                        netcdf_compression=setu['netcdf_compression'],
                                        netcdf_compression_level=setu['netcdf_compression_level'],
                                        netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
                    update_projection = False
                    new = False
                    if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data_full.shape))

            if verbosity > 1:
                print('Conversion took {:.1f} seconds'.format(time.time()-t0))
                print('Created {}'.format(ofile))

            if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
