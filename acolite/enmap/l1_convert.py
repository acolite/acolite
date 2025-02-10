## def l1_convert
## converts EnMAP file to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2022-09-20
## modifications: 2022-09-21 (QV) added band outputs
##                2022-11-16 (QV) added F0 reference
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-02-22 (QV) added computation of view angles
##                2024-04-16 (QV) use new gem NetCDF handling
##                2025-01-30 (QV) moved polygon limit
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use

def l1_convert(inputfile, output = None, settings = None):
    import numpy as np
    import scipy
    import datetime, dateutil.parser, os, copy
    import acolite as ac

    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    sensor = 'ENMAP_HSI'

    ## get sensor specific defaults
    setd = ac.acolite.settings.parse(sensor)
    ## set sensor default if user has not specified the setting
    for k in setd:
        if k not in ac.settings['user']: setu[k] = setd[k]
    ## end set sensor specific defaults

    verbosity = setu['verbosity']
    if output is None: output = setu['output']

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

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
        if (setu['limit'] is not None) and (dct is not None):
            dct_sub = ac.shared.projection_sub(dct, setu['limit'])
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

        ## set up oname (without directory or file type) and ofile (with directory and file type)
        oname  = '{}_{}'.format(gatts['sensor'],  time.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname,  gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## compute viewing angles
        ## 2024-02-22 QV
        roll  = -metadata['acrossOffNadirAngle']['center']
        pitch = -metadata['alongOffNadirAngle']['center']
        yaw   = -metadata['sceneAzimuthAngle']['center']
        roll  = roll / np.cos(np.radians(pitch))
        rot = scipy.spatial.transform.Rotation.from_euler('yxz', [pitch, roll, yaw], degrees=True)
        v_ = rot.apply(np.array([0, 0, 1]))
        cos = np.min((v_.dot(np.array([0, 0, 1])), 1.0))
        vza = np.degrees(np.arccos(cos))
        vaa = np.degrees(np.arctan2(-v_[0], -v_[1])) - 90
        if vaa < 0: vaa += 360
        gatts['vaa'] = vaa
        gatts['vza'] = vza

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
            if setu['polygon_clip']:
                clip_mask = ac.shared.polygon_crop(dct_prj, setu['polygon'], return_sub=False)
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

        ## set up output file
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gatts[k] for k in gatts}
        if dct_prj is not None:
            print('Computing and writing lat/lon')
            ## offset half pixels to compute center pixel lat/lon
            dct_prj['xrange'] = dct_prj['xrange'][0]+dct_prj['pixel_size'][0]/2, dct_prj['xrange'][1]-dct_prj['pixel_size'][0]/2
            dct_prj['yrange'] = dct_prj['yrange'][0]+dct_prj['pixel_size'][1]/2, dct_prj['yrange'][1]-dct_prj['pixel_size'][1]/2
            ## compute lat/lon
            lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel = False)
            print(lat.shape)
            gemo.write('lon', lon)
            lon = None
            gemo.write('lat', lat)
            lon = None

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
            print('Computing rhot_{} for {}'.format(bands[b]['wave_name'], gatts['oname']))
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

            if (setu['polygon_clip']): cdata_radiance[clip_mask] = np.nan

            if setu['output_lt']:
                ## write toa radiance
                gemo.write('Lt_{}'.format(bands[b]['wave_name']), cdata_radiance, ds_att = ds_att)

            ## compute reflectance
            cdata = cdata_radiance * (np.pi * gatts['se_distance'] * gatts['se_distance']) / (bands[b]['f0'] * mu0)
            cdata_radiance = None
            gemo.write('rhot_{}'.format(bands[b]['wave_name']), cdata, ds_att = ds_att)
            cdata = None
        cube, cube_swir = None, None
        gemo.close()
        ofiles.append(ofile)
    return(ofiles, setu)
