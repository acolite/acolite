## def l1_convert
## converts Worldview data to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-02-25
## modifications:

def l1_convert(inputfile,
               inputfile_swir = None,
               output = None,
               limit = None, sub = None,
               poly = None,
               output_geolocation = True,
               percentiles_compute = True,
               percentiles = (0,1,5,10,25,50,75,90,95,99,100),
               verbosity = 0, vname = ''):

    import os, glob, dateutil.parser, datetime, time
    import numpy as np
    import scipy.interpolate
    import acolite as ac

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)

    ## parse swir inputfile
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

    ## check if ROI polygon is given
    clip, clip_mask = False, None
    if poly is not None:
        if os.path.exists(poly):
            try:
                limit = ac.shared.polygon_limit(poly)
                print('Using limit from polygon envelope: {}'.format(limit))
                print('Not yet implemented for WorldView')
                clip = True
            except:
                print('Failed to import polygon {}'.format(poly))

    new = True
    ofile = None
    ofiles = []

    for fi, bundle in enumerate(inputfile):
        t0 = time.time()
        swir_bundle = inputfile_swir[fi] if inputfile_swir is not None else None

        ## parse the metadata
        if verbosity > 1: print('Importing metadata from {}'.format(bundle))
        metafile = glob.glob('{}/{}'.format(bundle,'*.XML'))
        if len(metafile)>0:
            metafile = metafile[0]
            meta = ac.worldview.metadata_parse(metafile)
        else:
            print('No metadata found for {}'.format(bundle))
            continue

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

        ## get F0 - not stricty necessary if using USGS reflectance
        f0 = ac.shared.f0_get()
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data'])*10, rsr)

        gatts = {'sensor':meta['sensor'], 'satellite':meta['satellite'],
                     'isodate':isodate, #'global_dims':global_dims,
                     'sza':sza, 'vza':vza, 'raa':raa, 'se_distance': se_distance,
                     'mus': np.cos(sza*(np.pi/180.)), 'acolite_file_type': 'L1R'}

        stime = dateutil.parser.parse(gatts['isodate'])
        oname = '{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'))
        if vname != '': oname+='_{}'.format(vname)

        ofile = '{}/{}_L1R.nc'.format(output, oname)
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
        ## get list of tiles in this bundle
        ntiles = len(meta['TILE_INFO'])
        for ti, tile_mdata in enumerate(meta['TILE_INFO']):
            tile = tile_mdata['FILENAME'].split('_')[1].split('-')[0]
            file = '{}/{}'.format(bundle,tile_mdata['FILENAME'])
            ## check if the files were named .TIF instead of .TIFF
            if not os.path.exists(file): file = file.replace('.TIFF', '.TIF')
            if not os.path.exists(file): continue
            try:
                ## get projection info from current tile
                dct_vnir = ac.shared.projection_read(file)
                ## set up dict
                if dct is None:
                    dct = {k: dct_vnir[k] for k in dct_vnir}
                else:
                    ## compute new ranges
                    dct['xrange'] = min(dct['xrange'][0], dct_vnir['xrange'][0]),\
                                    max(dct['xrange'][1], dct_vnir['xrange'][1])
                    dct['yrange'] = max(dct['yrange'][0], dct_vnir['yrange'][0]),\
                                    min(dct['yrange'][1], dct_vnir['yrange'][1])
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

        ## final scene dimensions
        if dct is not None:
            ## compute dimensions
            dct['xdim'] = int(np.round((dct['xrange'][1]-dct['xrange'][0]) / dct['pixel_size'][0]))
            dct['ydim'] = int(np.round((dct['yrange'][1]-dct['yrange'][0]) / dct['pixel_size'][1]))
            ## these should match the global dims from metadata
            if (global_dims[0] != dct['ydim']) |  (global_dims[1] != dct['xdim']):
                print('Global dims and projection size do not match')
                print(global_dims[1], dct['xdim'])
                print(global_dims[0], dct['ydim'])
            ## add projection to gatts
            for k in ['xrange', 'yrange', 'proj4_string', 'pixel_size', 'dimensions', 'zone']:
                if k in dct: gatts[k] = dct[k]

        ## write results to output file
        new=True

        ## write lat/lon
        if output_geolocation:
            if verbosity > 1: print('{} - Writing lat/lon'.format(datetime.datetime.now().isoformat()[0:19]))
            if dct is not None: ## compute from projection info
                lon, lat = ac.shared.projection_geo(dct, add_half_pixel=False)
                ac.output.nc_write(ofile, 'lat', lat, global_dims=global_dims, new=new, attributes=gatts)
                lat = None
                ac.output.nc_write(ofile, 'lon', lon)
                lon = None
            else: ## compute from corners given in metadata
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
                zlon = scipy.interpolate.interp2d(pcol, prow, plon, kind='linear')
                zlat = scipy.interpolate.interp2d(pcol, prow, plat, kind='linear')
                x = np.arange(1, 1+global_dims[1], 1)
                y = np.arange(1, 1+global_dims[0], 1)
                ac.output.nc_write(ofile, 'lat', zlat(x, y), global_dims=global_dims, new=new, attributes=gatts)
                ac.output.nc_write(ofile, 'lon', zlon(x, y))
                x = None
                y = None
                zlat = None
                zlon = None
            new = False
        ## end write lat/lon

        ## run through bands
        for b,band in enumerate(band_names):
            ## run through tiles in this bundle
            ntiles = len(meta['TILE_INFO'])
            for ti, tile_mdata in enumerate(meta['TILE_INFO']):
                tile = tile_mdata['FILENAME'].split('_')[1].split('-')[0]

                ## get tile offset
                offset = [int(tile_mdata['ULCOLOFFSET']), int(tile_mdata['ULROWOFFSET'])]
                if verbosity > 1: print('{} - Band {} Processing tile {}/{}'.format(datetime.datetime.now().isoformat()[0:19], band, ti+1, ntiles), tile, offset)

                file = '{}/{}'.format(bundle,tile_mdata['FILENAME'])
                ## check if the files were named .TIF instead of .TIFF
                if not os.path.exists(file): file = file.replace('.TIFF', '.TIF')
                if not os.path.exists(file): continue

                ## get tile from wv3 swir bundle if provided
                swir_file=None
                if swir_bundle is not None:
                    for tile_mdata_swir in swir_meta['TILE_INFO']:
                        if tile in tile_mdata_swir['FILENAME']:
                            swir_file = '{}/{}'.format(swir_bundle,tile_mdata_swir['FILENAME'])
                        if not os.path.exists(swir_file):
                            continue
                ## end swir bundle

                if 'SWIR' not in band:
                    bt = [bt for bt in meta['BAND_INFO'] if meta['BAND_INFO'][bt]['name'] == band][0]
                    bd = meta['BAND_INFO'][bt]
                    d = ac.shared.read_band(file, idx=bd['index'], sub=sub)
                else:
                    if swir_file is None:
                        swir_file='{}'.format(file)
                        swir_meta = meta.copy()
                    bt = [bt for bt in swir_meta['BAND_INFO'] if swir_meta['BAND_INFO'][bt]['name'] == band][0]
                    bd = swir_meta['BAND_INFO'][bt]
                    d = ac.shared.read_band(swir_file, idx=bd['index'], sub=sub)

                ## track mask
                nodata = d == np.uint16(0)
                ## convert to float and scale to TOA reflectance
                d = d.astype(np.float32)
                d *= float(meta['BAND_INFO'][bt]['ABSCALFACTOR'])/float(meta['BAND_INFO'][bt]['EFFECTIVEBANDWIDTH'])
                d *= (np.pi * gatts['se_distance']**2) / (f0_b[band]/10. * gatts['mus'])
                ## apply mask
                d[nodata] = np.nan

                ## make new data full array
                if ti == 0: data_full = np.zeros(global_dims) + np.nan

                ## add in data
                data_full[offset[1]:offset[1]+d.shape[0], offset[0]:offset[0]+d.shape[1]] = d
                d = None

            ## set up dataset attributes
            ds = 'rhot_{}'.format(waves_names[band])
            ds_att = {'wavelength':waves_mu[band]*1000}
            if percentiles_compute:
                ds_att['percentiles'] = percentiles
                ds_att['percentiles_data'] = np.nanpercentile(data_full, percentiles)

            ## write to netcdf file
            if verbosity > 1: print('{} - Converting bands: Writing {} ({})'.format(datetime.datetime.now().isoformat()[0:19], ds, data_full.shape))
            ac.output.nc_write(ofile, ds, data_full, attributes = gatts, new = new, dataset_attributes = ds_att)
            if verbosity > 1: print('{} - Converting bands: Wrote {} ({})'.format(datetime.datetime.now().isoformat()[0:19], ds, data_full.shape))
            new = False
            data_full = None

        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles)
