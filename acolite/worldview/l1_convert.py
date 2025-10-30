## def l1_convert
## converts Worldview data to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-02-25
## modifications: 2021-12-31 (QV) new handling of settings
##                2022-01-04 (QV) added netcdf compression
##                2022-11-14 (QV) added subsetting of projected data
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-17 (QV) use new gem NetCDF handling
##                2025-01-28 (QV) switch to LinearNDInterpolator, added meshgrid
##                2025-01-30 (QV) moved polygon limit
##                2025-02-02 (QV) removed percentiles, switched to ac.shared.interp2d
##                2025-02-03 (QV) updated tile merging for projected data and user set region of interest
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use
##                2025-03-01 (QV) added bundle_test to get metafile, added PGC identification
##                                replaced convert_atmospherically_corrected by worldview_convert_l2 setting
##                2025-04-07 (QV) change worldview_convert_l2 to convert_l2
##                2025-09-15 (QV) added rpc_use = False

def l1_convert(inputfile, output = None, inputfile_swir = None, settings = None):

    import os, glob, dateutil.parser, datetime, time
    import numpy as np
    import scipy.interpolate
    import acolite as ac

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

    ## parse swir inputfile
    if setu['inputfile_swir'] is not None: inputfile_swir = setu['inputfile_swir']
    if inputfile_swir is not None:
        if type(inputfile_swir) != list:
            if type(inputfile_swir) == str:
                inputfile_swir = inputfile_swir.split(',')
            else:
                inputfile_swir = list(inputfile_swir)
        nscenes_swir = len(inputfile_swir)
        if nscenes_swir != nscenes:
            print('Different number of scenes and SWIR scenes given.')
            return()
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    new = True
    ofile = None
    ofiles = []

    for fi, bundle in enumerate(inputfile):
        sub = None
        t0 = time.time()
        swir_bundle = inputfile_swir[fi] if inputfile_swir is not None else None

        ## parse the metadata
        if verbosity > 1: print('Importing metadata from {}'.format(bundle))
        metafile = ac.worldview.bundle_test(bundle)
        if metafile is None:
            print('No metadata found for {}'.format(bundle))
            continue
        else:
            meta = ac.worldview.metadata_parse(metafile)
        sensor = meta['sensor']

        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults

        verbosity = setu['verbosity']
        if output is None: output = setu['output']

        ## test if we need to do an atmospheric correction
        atmospherically_corrected = False
        if ('RADIOMETRICLEVEL' in meta) & ('RADIOMETRICENHANCEMENT' in meta):
            if meta['RADIOMETRICENHANCEMENT'] in ['ACOMP']:
                print('Image {} is already corrected by supplier.'.format(bundle))
                print('RADIOMETRICLEVEL: {} RADIOMETRICENHANCEMENT: {}'.format(meta['RADIOMETRICLEVEL'], meta['RADIOMETRICENHANCEMENT']))
                atmospherically_corrected = True
                if not setu['convert_l2']: continue

        ## test if we have PGC bundle
        pgc_bundle = False
        if meta['PGC']:
            pgc_bundle = True
            if meta['PGC_STRETCH'] in ['mr']:
                print('Image {} is already corrected by supplier.'.format(bundle))
                print('PGC_STRETCH: {}'.format(meta['PGC_STRETCH']))
                atmospherically_corrected = True
                if not setu['convert_l2']: continue

        ## parse the metadata
        if swir_bundle is not None:
            if verbosity > 1: print('Importing metadata from {}'.format(swir_bundle))
            swir_metafile = glob.glob('{}/{}'.format(swir_bundle,'*.XML'))
            if len(swir_metafile)>0:
                swir_metafile = swir_metafile[0]
                swir_meta = ac.worldview.metadata_parse(swir_metafile)
            else:
                print('No metadata found for {}'.format(swir_bundle))
                continue

        band_names = [meta['BAND_INFO'][b]['name'] for b in list(meta['BAND_INFO'].keys())]
        if len(band_names) == 0:
            print('Could not identify required bands from {}'.format(bundle))
            print('Is this a multispectral WorldView image?')
            continue
        if swir_bundle is not None: band_names += [swir_meta['BAND_INFO'][b]['name'] for b in list(swir_meta['BAND_INFO'].keys())]

        ## get observation geometry
        raa = abs(float(meta['MEANSUNAZ']) - float(meta['MEANSATAZ']))
        while raa >= 180.: raa = np.abs(raa-360)
        sza = 90. - float(meta['MEANSUNEL'])
        vza = 90. - float(meta['MEANSATEL'])

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

        gains = None
        if setu['gains']:
            if (len(setu['gains_toa']) == len(rsr_bands)) &\
               (len(setu['offsets_toa']) == len(rsr_bands)):
               gains = {}
               for bi, band in enumerate(rsr_bands):
                   gains[band] = {'gain': float(setu['gains_toa'][bi]),
                                'offset': float(setu['offsets_toa'][bi])}
            else:
                print('Use of gains requested, but provided number of gain ({}) or offset ({}) values does not match number of bands in RSR ({})'.format(len(setu['gains_toa']), len(setu['offsets_toa']), len(rsr_bands)))
                print('Provide gains in band order: {}'.format(','.join(rsr_bands)))

        ## get F0 - not stricty necessary if using USGS reflectance
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data'])*10, rsr)

        gatts = {'sensor':meta['sensor'], 'satellite':meta['satellite'],
                     'isodate':isodate, #'global_dims':global_dims,
                     'sza':sza, 'vza':vza, 'raa':raa, 'se_distance': se_distance,
                     'mus': np.cos(sza*(np.pi/180.)), 'acolite_file_type': 'L1R'}

        if atmospherically_corrected: gatts['acolite_file_type'] = 'converted'

        stime = dateutil.parser.parse(gatts['isodate'])

        ## set up oname (without directory or file type) and ofile (with directory and file type)
        oname = '{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## add band info to gatts
        for b in rsr_bands:
            gatts['{}_wave'.format(b)] = waves_mu[b]*1000
            gatts['{}_name'.format(b)] = waves_names[b]
            gatts['{}_f0'.format(b)] = f0_b[b]

        ## global scene dimensions from metadata
        global_dims = [int(meta['NUMROWS']),int(meta['NUMCOLUMNS'])]

        ## get projection info from tiles
        dct = None
        nc_projection = None
        warp_to = None
        rpc_use = False

        ## get list of tiles in this bundle
        ntiles = len(meta['TILE_INFO'])
        for ti, tile_mdata in enumerate(meta['TILE_INFO']):
            try:
                tile = tile_mdata['FILENAME'].split('_')[1].split('-')[0]
            except:
                tile = ''

            file = '{}/{}'.format(bundle,tile_mdata['FILENAME'])

            ## check if the files were named .TIF instead of .TIFF
            if not os.path.exists(file): file = file.replace('.TIFF', '.TIF')
            if not os.path.exists(file): continue

            ## get projection info from current tile
            try:
                dct_vnir = ac.shared.projection_read(file)
                if setu['verbosity'] > 5:
                    print('Read projection info from {}'.format(file))
                    print(dct_vnir)
            except:
                dct_vnir = None
                rpc_use = True

            # ## shall we reproject the inputfile?
            # ## to add swir bundle
            # reproject = setu['reproject_inputfile_force'] | (setu['reproject_inputfile'] & (dct_vnir == None))
            # if reproject:
            #     setu['export_geotiff_match_file'] = None ## OG extent will no longer match
            #     bn = os.path.basename(file)
            #     bn, ext = os.path.splitext(bn)
            #     ifile = '{}'.format(file)
            #     rfile = '{}/{}_reprojected{}'.format(output, bn, ext)
            #     print(ifile)
            #     print(rfile)
            #     #from osgeo import gdal
            #     #print('Reprojecting {} to {}'.format(ifile, file))
            #     #dso = gdal.Warp(file, ifile) ## warp to gdal defaults
            #     file, (dimxo, dimyo) = ac.shared.warp_inputfile(ifile, target=rfile)
            #     ## update global dims
            #     global_dims = dimyo, dimxo
            #     print(file)
            #     #dso = None
            #     dct_vnir = ac.shared.projection_read(file)

            ## get projection info
            try:
                ## set up dict
                if dct is None:
                    dct = {k: dct_vnir[k] for k in dct_vnir}
                else:
                    ## compute new ranges
                    dct['xrange'] = min(dct['xrange'][0], dct_vnir['xrange'][0]),\
                                    max(dct['xrange'][1], dct_vnir['xrange'][1])
                    dct['yrange'] = max(dct['yrange'][0], dct_vnir['yrange'][0]),\
                                    min(dct['yrange'][1], dct_vnir['yrange'][1])
                    ## update dimensions
                    dct['xdim'] = np.round((dct['xrange'][1] - dct['xrange'][0]) / dct['pixel_size'][0]).astype(int)
                    dct['ydim'] = np.round((dct['yrange'][1] - dct['yrange'][0]) / dct['pixel_size'][1]).astype(int)
                    dct['dimensions'] = (dct['xdim'], dct['ydim'])
                    if setu['verbosity'] > 5:
                        print('Updated projection info from {}'.format(file))
                        print(dct)

                ## get projection info from swir file - not used atm, check if SWIR band projection matches?
                #swir_file=None
                #if swir_bundle is not None:
                #    for tile_mdata_swir in swir_meta['TILE_INFO']:
                #        if tile in tile_mdata_swir['FILENAME']:
                #            swir_file = '{}/{}'.format(swir_bundle,tile_mdata_swir['FILENAME'])
                #        if not os.path.exists(swir_file): continue
                #if swir_file is not None: dct_swir = ac.shared.projection_read(swir_file)
            except:
                if verbosity > 1: print('Could not determine projection from {}'.format(file))
                pass

        ## set up limit and projection dct so we can warp to target projection when reading data
        if (dct is None) & (setu['worldview_reproject']):
            if setu['limit'] is not None:
                lim = setu['limit']
            else:
                bt = list(meta['BAND_INFO'].keys())[0]
                lons = [meta['BAND_INFO'][bt][k] for k in meta['BAND_INFO'][bt] if 'LON' in k]
                lats = [meta['BAND_INFO'][bt][k] for k in meta['BAND_INFO'][bt] if 'LAT' in k]
                lim = [min(lats), min(lons), max(lats), max(lons)]
            dct, nc_projection, warp_to = ac.shared.projection_setup(lim, setu['worldview_reproject_resolution'], \
                                                                          res_method=setu['worldview_reproject_method'])
            global_dims = dct['ydim'], dct['xdim']

        ## final scene dimensions
        if dct is not None:
            ## if we have dct and limit we can subset
            if setu['limit'] is not None:
                dct_sub = ac.shared.projection_sub(dct, setu['limit'], four_corners=True)
                if dct_sub['out_lon']:
                    if verbosity > 1: print('Longitude limits outside {}'.format(bundle))
                    continue
                if dct_sub['out_lat']:
                    if verbosity > 1: print('Latitude limits outside {}'.format(bundle))
                    continue
                sub = dct_sub['sub']
                dct = {k:dct_sub[k] for k in dct_sub}
                global_dims = dct['ydim'], dct['xdim']

            ## create warp_to dataset
            warp_to = ac.shared.projection_warp_to(dct, res_method = setu['warp_resampling_method'])

            ## compute dimensions
            print(dct['xdim'], dct['ydim'])
            dct['xdim'] = int(np.round((dct['xrange'][1]-dct['xrange'][0]) / dct['pixel_size'][0]))
            dct['ydim'] = int(np.round((dct['yrange'][1]-dct['yrange'][0]) / dct['pixel_size'][1]))
            print(dct['xdim'], dct['ydim'])

            print(dct['xrange'], dct['yrange'])

            ## these should match the global dims from metadata
            if setu['limit'] is None:
                if (global_dims[0] != dct['ydim']) |  (global_dims[1] != dct['xdim']):
                    print('Global dims and projection size do not match')
                    print(global_dims[1], dct['xdim'])
                    print(global_dims[0], dct['ydim'])
            ## add projection to gatts
            for k in ['xrange', 'yrange', 'proj4_string', 'pixel_size', 'dimensions', 'zone']:
                if k in dct: gatts[k] = dct[k]

            if nc_projection is None: nc_projection = ac.shared.projection_netcdf(dct, add_half_pixel = False)
            gatts['projection_key'] = [k for k in nc_projection if k not in ['x', 'y']][0]

        ## write results to output file
        new = True
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gatts[k] for k in gatts}
        gemo.nc_projection = nc_projection
        new = False

        ## write lat/lon
        if setu['output_geolocation']:
            if verbosity > 1: print('{} - Writing lat/lon'.format(datetime.datetime.now().isoformat()[0:19]))
            if dct is not None: ## compute from projection info
                print('lat/lon computed from projection info')
                lon, lat = ac.shared.projection_geo(dct, add_half_pixel=False)
                gemo.write('lon', lon)
                lon = None
                gemo.write('lat', lat)
                lat = None
            else: ## compute from corners given in metadata
                print('lat/lon interpolated from metadata corners')
                pcol = [0, global_dims[1], global_dims[1], 0]
                prow = [0, 0, global_dims[0], global_dims[0]]
                plon = []
                plat = []
                band_tag = list(meta['BAND_INFO'].keys())[0]
                for bk in ['UL', 'UR', 'LR', 'LL']:
                        k = '{}{}'.format(bk, 'LON')
                        plon.append(meta['BAND_INFO'][band_tag][k])
                        k = '{}{}'.format(bk, 'LAT')
                        plat.append(meta['BAND_INFO'][band_tag][k])

                ## set up interpolator
                zlon = ac.shared.interp2d((pcol, prow), plon)
                zlat = ac.shared.interp2d((pcol, prow), plat)
                x = np.arange(1, 1+global_dims[1], 1)
                y = np.arange(1, 1+global_dims[0], 1)
                X, Y = np.meshgrid(x, y)
                gemo.write('lon', zlon(X, Y))
                gemo.write('lat', zlat(X, Y))
                del X, Y
                del x, y, zlat, zlon
        ## end write lat/lon

        ## run through bands
        for b,band in enumerate(band_names):
            data_full = None

            ## run through tiles in this bundle
            ntiles = len(meta['TILE_INFO'])
            for ti, tile_mdata in enumerate(meta['TILE_INFO']):
                try:
                    tile = tile_mdata['FILENAME'].split('_')[1].split('-')[0]
                except:
                    tile = ''

                ## continue if subsetting is not correct and warp_to is not present
                if (warp_to is None) & (sub is not None):
                    if sub[1] <= sub[0]: continue
                    if sub[3] <= sub[2]: continue

                ## get tile offset
                offset = [int(tile_mdata['ULCOLOFFSET']), int(tile_mdata['ULROWOFFSET'])]
                if verbosity > 1: print('{} - Band {} Processing tile {}/{}'.format(datetime.datetime.now().isoformat()[0:19], band, ti+1, ntiles), tile, offset)

                file = '{}/{}'.format(bundle,tile_mdata['FILENAME'])
                ## check if the files were named .TIF instead of .TIFF
                if not os.path.exists(file): file = file.replace('.TIFF', '.TIF')
                if not os.path.exists(file): continue

                # ## to add swir bundle?
                # if reproject:
                #     bn = os.path.basename(file)
                #     bn, ext = os.path.splitext(bn)
                #     ifile = '{}'.format(file)
                #     file = '{}/{}_reprojected{}'.format(output, bn, ext)
                #     print(file)

                ## get tile from wv3 swir bundle if provided
                swir_file=None
                if swir_bundle is not None:
                    for tile_mdata_swir in swir_meta['TILE_INFO']:
                        if tile in tile_mdata_swir['FILENAME']:
                            swir_file = '{}/{}'.format(swir_bundle,tile_mdata_swir['FILENAME'])
                        if not os.path.exists(swir_file):
                            continue
                ## end swir bundle

                ## get band scaling factors cf
                if 'SWIR' not in band:
                    bt = [bt for bt in meta['BAND_INFO'] if meta['BAND_INFO'][bt]['name'] == band][0]
                    d = ac.shared.read_band(file, idx=meta['BAND_INFO'][bt]['index'], sub = sub, warp_to = warp_to, rpc_use = rpc_use)
                    cf = float(meta['BAND_INFO'][bt]['ABSCALFACTOR'])/float(meta['BAND_INFO'][bt]['EFFECTIVEBANDWIDTH'])
                else:
                    if swir_file is None:
                        swir_file='{}'.format(file)
                        swir_meta = meta.copy()
                    bt = [bt for bt in swir_meta['BAND_INFO'] if swir_meta['BAND_INFO'][bt]['name'] == band][0]
                    d = ac.shared.read_band(swir_file, idx=swir_meta['BAND_INFO'][bt]['index'], sub = sub, warp_to = warp_to, rpc_use = rpc_use)
                    cf = float(swir_meta['BAND_INFO'][bt]['ABSCALFACTOR'])/float(swir_meta['BAND_INFO'][bt]['EFFECTIVEBANDWIDTH'])

                ## skip if one dimension is 0
                if (d.shape[0] == 0) | (d.shape[1] == 0):
                    continue

                if cf <=0:
                    print('Warning DN scaling factor is <0, this will give bad TOA radiances/reflectances.')
                    if 'RADIOMETRICENHANCEMENT' in meta:
                        print('Data has been enhanced by the provider: {}'.format(meta['RADIOMETRICENHANCEMENT']))

                ## track mask
                dtype = d.dtype
                if dtype == np.dtype('uint8'):
                    nodata = d == np.uint8(0)
                elif dtype == np.dtype('uint16'):
                    nodata = d == np.uint16(0)
                elif dtype == np.dtype('float32'):
                    nodata = d == np.float32(0)

                ## override cf for PGC bundle
                if pgc_bundle:
                    if dtype == np.dtype('uint8'): cf = 1/200.
                    elif dtype == np.dtype('uint16'): cf = 1/2000.
                    else : cf = 1.0
                    print('PGC scaling factor for stretch {}, dtype {}: {:.4f}'.format(meta['PGC_STRETCH'], dtype, cf))
                    cf *= (gatts['se_distance']**2) / gatts['mus']
                    print('PGC scaling factor with sun earth distance and zenith angle: {:.4f}'.format(cf))

                ## convert to float and scale to TOA reflectance
                d = d.astype(np.float32) * cf

                ## convert from original MAXAR bundle
                if not pgc_bundle:
                    if (gains != None) & (setu['gains_parameter'] == 'radiance'):
                        print('Applying gain {} and offset {} to TOA radiance for band {}'.format(gains[band]['gain'], gains[band]['offset'], band))
                        d = gains[band]['gain'] * d + gains[band]['offset']
                    d *= (np.pi * gatts['se_distance']**2) / (f0_b[band]/10. * gatts['mus'])
                    if (gains != None) & (setu['gains_parameter'] == 'reflectance'):
                        print('Applying gain {} and offset {} to TOA reflectance for band {}'.format(gains[band]['gain'], gains[band]['offset'], band))
                        d = gains[band]['gain'] * d + gains[band]['offset']

                ## apply mask
                d[nodata] = np.nan

                ## make new data full array for the current band
                if data_full is None: data_full = np.zeros(global_dims) + np.nan

                ## add in data using offset if shape does not match
                if d.shape != data_full.shape:
                    data_full[offset[1]:offset[1]+d.shape[0], offset[0]:offset[0]+d.shape[1]] = d
                else:
                    dsub = np.where(np.isnan(data_full))
                    data_full[dsub] = d[dsub]
                d = None
            if data_full is None: continue

            ## set up dataset attributes
            ds = 'rhot_{}'.format(waves_names[band])
            if atmospherically_corrected: ds = ds.replace('rhot_', 'rhos_acomp_')

            ds_att = {'wavelength': waves_mu[band]*1000, 'band_name': band, 'f0': f0_b[band]/10.}
            if gains != None:
                ds_att['gain'] = gains[band]['gain']
                ds_att['offset'] = gains[band]['offset']
                ds_att['gains_parameter'] = setu['gains_parameter']

            ## write to netcdf file
            if verbosity > 1: print('{} - Converting bands: Writing {} ({})'.format(datetime.datetime.now().isoformat()[0:19], ds, data_full.shape))
            gemo.write(ds, data_full, ds_att = ds_att)
            if verbosity > 1: print('{} - Converting bands: Wrote {} ({})'.format(datetime.datetime.now().isoformat()[0:19], ds, data_full.shape))
            data_full = None

        gemo.close()
        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
