## def l1_convert
## converts EnMAP file to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2022-09-20
## modifications: 2022-09-21 (QV) added band outputs
##                2022-11-16 (QV) added F0 reference

def l1_convert(inputfile, output = None, settings = {}, verbosity = 5):
    import numpy as np
    import datetime, dateutil.parser, os, copy
    import acolite as ac

    ## parse sensor specific settings
    setu = ac.acolite.settings.parse('ENMAP_HSI', settings=settings)
    vname = setu['region_name']
    output_lt = setu['output_lt']
    if output is None: output = setu['output']
    verbosity = setu['verbosity']
    poly = setu['polygon']
    limit = setu['limit']

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

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

    ## get F0 for radiance -> reflectance computation
    f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])

    ofiles = []
    for bundle in inputfile:
        bundle_files = ac.enmap.bundle_test(bundle)
        metadata, band_data = ac.enmap.metadata_parse(bundle_files['METADATA'])

        satellite = metadata['mission'].upper()
        sensor = metadata['sensor'].upper()

        imagefile, imagefile_swir = None, None
        if 'SPECTRAL_IMAGE' in bundle_files:
            imagefile = bundle_files['SPECTRAL_IMAGE']
        elif 'SPECTRAL_IMAGE_VNIR' in bundle_files:
            imagefile = bundle_files['SPECTRAL_IMAGE_VNIR']
        if 'SPECTRAL_IMAGE_SWIR' in bundle_files:
            imagefile_swir = bundle_files['SPECTRAL_IMAGE_SWIR']

        ## set up projection
        warp_to, dct_prj, sub = None, None, None
        try:
            ## get projection from image
            dct = ac.shared.projection_read(imagefile)
        except:
            print('Could not determine image projection')
            dct = None

        ## find crop
        if (limit is not None) and (dct is not None):
            dct_sub = ac.shared.projection_sub(dct, limit)
            if dct_sub['out_lon']:
                if verbosity > 1: print('Longitude limits outside {}'.format(bundle))
                continue
            if dct_sub['out_lat']:
                if verbosity > 1: print('Latitude limits outside {}'.format(bundle))
                continue
            sub = dct_sub['sub']

        if dct is not None:
            if sub is None:
                dct_prj = {k:dct[k] for k in dct}
            else:
                dct_prj = {k:dct_sub[k] for k in dct_sub}

                ## updated 2022-03-28
                xyr = [min(dct_prj['xrange']),
                       min(dct_prj['yrange']),
                       max(dct_prj['xrange']),
                       max(dct_prj['yrange']),
                       dct_prj['proj4_string']]

                ## warp settings for read_band
                res_method = 'near'
                warp_to = (dct_prj['proj4_string'], xyr, dct_prj['pixel_size'][0],dct_prj['pixel_size'][1], res_method)


        ## date and time
        stime = dateutil.parser.parse(metadata['startTime'])
        etime = dateutil.parser.parse(metadata['stopTime'])
        otime = (etime-stime).seconds
        time = stime + datetime.timedelta(seconds=otime/2)
        doy = time.strftime('%j')
        se_distance = ac.shared.distance_se(doy)


        ## collect global attributes
        gatts = {}
        gatts['acolite_file_type'] = 'L1R'
        gatts['isodate'] = time.isoformat()
        gatts['sensor'] = '{}_{}'.format(metadata['mission'], metadata['sensor']).upper()
        gatts['doy'] = doy
        gatts['se_distance'] = se_distance
        obase  = '{}_{}_L1R'.format(gatts['sensor'],  time.strftime('%Y_%m_%d_%H_%M_%S'))
        gatts['obase'] = obase

        ## compute viewing azimuth
        ## uncertain, check with EnMAP team
        #orientation_angle = metadata['sceneAzimuthAngle']['center']
        #ang_across = metadata['acrossOffNadirAngle']['center']
        #ang_along = metadata['alongOffNadirAngle']['center']
        #vaa = np.mod(orientation_angle - np.degrees(np.arctan2(np.tan(np.radians(ang_across)),
        #                                                       np.tan(np.radians(ang_along)))),360)
        #gatts['vaa'] = vaa
        #gatts['vza'] = 0

        ## edited 2022-09-22
        gatts['vaa'] = metadata['sceneAzimuthAngle']['center']
        gatts['vza'] = np.abs(metadata['acrossOffNadirAngle']['center'])

        gatts['saa'] = metadata['sunAzimuthAngle']['center']
        gatts['sza'] = 90-metadata['sunElevationAngle']['center']

        gatts['raa'] = np.abs(gatts['saa']-gatts['vaa'])
        while gatts['raa'] > 180: gatts['raa'] = abs(gatts['raa']-360)

        mu0 = np.cos(gatts['sza']*(np.pi/180))
        muv = np.cos(gatts['vza']*(np.pi/180))

        ## band information
        gatts['band_waves'] = [band_data[b]['wavelengthCenterOfBand'] for b in band_data]
        gatts['band_widths'] = [band_data[b]['FWHMOfBand'] for b in band_data]

        ## add projection info
        if dct_prj is not None:
            pkeys = ['xrange', 'yrange', 'proj4_string', 'pixel_size', 'zone']
            for k in pkeys:
                if k in dct_prj: gatts[k] = copy.copy(dct_prj[k])

            ## if we are clipping to a given polygon get the clip_mask here
            if clip:
                clip_mask = ac.shared.polygon_crop(dct_prj, poly, return_sub=False)
                clip_mask = clip_mask.astype(bool) == False

        ## band rsr and f0
        rsr = ac.shared.rsr_hyper(gatts['band_waves'], gatts['band_widths'])
        f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], rsr)

        ## make bands dataset
        bands = {}
        for bi, b in enumerate(rsr):
            cwave = gatts['band_waves'][bi]
            swave = '{:.0f}'.format(cwave)
            bands[b]= {'wave':cwave, 'wavelength':cwave, 'wave_mu':cwave/1000.,
                       'wave_name':'{:.0f}'.format(cwave),
                       'width': gatts['band_widths'][bi],
                       'rsr': rsr[b],'f0': f0d[b]}

        if output is None:
            odir = os.path.dirname(imagefile)
        else:
            odir = output
        if not os.path.exists(odir): os.makedirs(odir)
        ofile = '{}/{}.nc'.format(odir, obase)

        new = True
        if dct_prj is not None:
            print('Computing and writing lat/lon')
            ## offset half pixels to compute center pixel lat/lon
            dct_prj['xrange'] = dct_prj['xrange'][0]+dct_prj['pixel_size'][0]/2, dct_prj['xrange'][1]-dct_prj['pixel_size'][0]/2
            dct_prj['yrange'] = dct_prj['yrange'][0]+dct_prj['pixel_size'][1]/2, dct_prj['yrange'][1]-dct_prj['pixel_size'][1]/2
            ## compute lat/lon
            lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel = False)
            print(lat.shape)
            ac.output.nc_write(ofile, 'lat', lat, new = new, attributes = gatts,
                                netcdf_compression=setu['netcdf_compression'],
                                netcdf_compression_level=setu['netcdf_compression_level'])
            lat = None
            ac.output.nc_write(ofile, 'lon', lon,
                                netcdf_compression=setu['netcdf_compression'],
                                netcdf_compression_level=setu['netcdf_compression_level'])
            lon = None
            new = False

        ## read data cube (faster)
        read_cube = True
        if read_cube:
            print('Reading EnMAP image cube')
            cube = ac.shared.read_band(imagefile, sub = sub, warp_to = warp_to).astype(np.float32)
            cube[cube==0.0] = np.nan
            print('Read EnMAP image cube {}'.format(cube.shape))
            if imagefile_swir != None:
                print('Reading EnMAP SWIR image cube')
                cube_swir = ac.shared.read_band(imagefile_swir, sub = sub, warp_to = warp_to).astype(np.float32)
                cube_swir[cube_swir==0.0] = np.nan
                print('Read EnMAP SWIR image cube {}'.format(cube_swir.shape))

        ## write TOA data
        for bi, b in enumerate(bands):
            bname = '{}'.format(bi+1)
            print('Computing rhot_{} for {}'.format(bands[b]['wave_name'], gatts['obase']))
            ds_att = {k: bands[b][k] for k in bands[b] if k not in ['rsr']}

            ## read data
            if read_cube:
                if (imagefile_swir != None) & (bi >= 91):
                    cdata_radiance = 1.0 * cube_swir[bi - 91, :, :]
                else:
                    cdata_radiance = 1.0 * cube[bi, :, :]
            else:
                if (imagefile_swir != None) & (bi >= 91):
                    cdata_radiance = ac.shared.read_band(imagefile_swir, bi-91+1, sub=sub, warp_to = warp_to).astype(np.float32)
                else:
                    cdata_radiance = ac.shared.read_band(imagefile, bi+1, sub=sub, warp_to = warp_to).astype(np.float32)
                cdata_radiance[cdata_radiance == 0] = np.nan

            ## compute radiance
            cdata_radiance *= band_data[bname]['GainOfBand']
            cdata_radiance += band_data[bname]['OffsetOfBand']
            cdata_radiance *= 1000 ## gains are W/m2/sr/nm so factor 1000 needed to get mW/m2/sr/nm

            if (clip) & (clip_mask is not None): cdata_radiance[clip_mask] = np.nan

            if output_lt:
                ## write toa radiance
                ac.output.nc_write(ofile, 'Lt_{}'.format(bands[b]['wave_name']), cdata_radiance,
                                            attributes = gatts, dataset_attributes = ds_att, new = new,
                                            netcdf_compression=setu['netcdf_compression'],
                                            netcdf_compression_level=setu['netcdf_compression_level'],
                                            netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
                new = False

            ## compute reflectance
            cdata = cdata_radiance * (np.pi * gatts['se_distance'] * gatts['se_distance']) / (bands[b]['f0'] * mu0)
            cdata_radiance = None

            ac.output.nc_write(ofile, 'rhot_{}'.format(bands[b]['wave_name']), cdata,\
                                            attributes = gatts, dataset_attributes = ds_att, new = new,
                                            netcdf_compression=setu['netcdf_compression'],
                                            netcdf_compression_level=setu['netcdf_compression_level'],
                                            netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
            cdata = None
            new = False
        cube, cube_swir = None, None
        ofiles.append(ofile)
    return(ofiles, setu)
