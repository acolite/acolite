## def l1_convert
## convert Huanjing image to L1R ACOLITE NetCDF file
##
## written by Quinten Vanhellemont, RBINS
## 2025-07-29
## modifications: 2025-07-30 (QV) added radiometric calibration, extent test

def l1_convert(inputfile, output = None, settings = None):
    import os, datetime, zoneinfo, dateutil.parser
    import numpy as np
    import acolite as ac

    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if setu['verbosity'] > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ## get radiometric calibration, converted from Excel Jianzhong 2025-07-30
    radiometric_calibration = ac.huanjing.get_radiometric_calibration()

    ## run through inputfiles
    ofiles = []
    for bundle in inputfile:
        imagefile, metafile = ac.huanjing.bundle_test(bundle)
        if (imagefile is None) | (metafile is None): continue

        ## parse metadata
        metadata = ac.huanjing.metadata_parse(metafile)
        sensor = '{}_{}'.format(metadata['SatelliteID'], metadata['SensorID'])

        rsrd = ac.shared.rsr_dict(sensor)[sensor]
        bands = rsrd['rsr_bands']

        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults
        if output is None: output = setu['output']

        resolution = setu['default_projection_resolution']

        ## get F0 for radiance -> reflectance computation
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], rsrd['rsr'])

        ## compute date time
        start_dt = dateutil.parser.parse(metadata['StartTime'])
        end_dt = dateutil.parser.parse(metadata['EndTime'])
        dt_ = start_dt + (end_dt - start_dt)
        dt_ = dt_.replace(tzinfo=zoneinfo.ZoneInfo('Asia/Shanghai')) ## China Standard Time?
        dt = dt_.astimezone(zoneinfo.ZoneInfo('UTC')) ## convert to UTC
        isodate = dt.isoformat()

        # ## find image subset
        sub, rpc_dem, warp_to = None, None, None
        dct, nc_projection = None, None

        ## get scene average geometry
        centre_lon = np.atleast_1d(float(metadata['CenterLongitude']))
        centre_lat = np.atleast_1d(float(metadata['TopLeftLatitude']))
        lat_range = [float(metadata[k]) for k in metadata if 'Latitude' in k]
        lon_range = [float(metadata[k]) for k in metadata if 'Longitude' in k]
        scene_limit = [min(lat_range), min(lon_range), max(lat_range), max(lon_range)]
        print(scene_limit)

        ## do reprojection if LEVEL1A
        if (metadata['ProductLevel'] == 'LEVEL1A'):
            if (setu['reproject_inputfile']):
                ## create extent of scene
                if setu['limit'] is None:
                    limit = [l for l in scene_limit]
                else:
                    limit = setu['limit']
                    if ((limit[0] > scene_limit[2]) |  (limit[2] < scene_limit[0]) | (limit[1] > scene_limit[3]) |  (limit[3] < scene_limit[1])):
                        print('Requested limit ({}) out of scene limit ({})'.format(limit, scene_limit))
                        continue

                print('Creating projection target for limit {}'.format(limit))
                dct, nc_projection, warp_to = ac.shared.projection_setup(limit, resolution, res_method = 'bilinear', utm = True, epsg = None, add_half_pixel=False)
                print(dct, warp_to)
                if setu['reproject_inputfile_dem']:
                    rpc_dem = ac.dem.copernicus.copernicus_dem_rpc(dct, output = output)
                    print(rpc_dem)
            elif (setu['limit'] is not None):
                print('Processing of a subset requested for LEVEL1A data without reprojection, which is not supported.')
                print('Processing of a subset for LEVEL1A requires reproject_inputfile=True')
                continue
        else:
            print('ProductLevel {} to be configured'.format(metadata['ProductLevel']))
            continue

        ## scene average geometry
        saa = float(metadata['SolarAzimuth'])
        sza = float(metadata['SolarZenith'])
        vaa = float(metadata['SatelliteAzimuth'])
        vza = float(metadata['SatelliteZenith'])

        ## relative azimuth
        raa = np.abs(saa-vaa)
        if raa>180: raa = 360 - raa
        print(saa, sza, vaa, vza, raa)

        ## compute sun geometry - currently only for sun earth distance
        print('Computing sun position for {}'.format(dt.isoformat()))
        spos = ac.shared.sun_position(isodate, centre_lon, centre_lat)
        se_distance = spos['distance']
        del spos

        ## cosine sun zenith angle
        mus = np.cos(np.radians(sza))

        ## set up gatts
        gatts = {}
        gatts['sensor'] = sensor
        gatts['isodate'] = dt.isoformat()
        gatts['acolite_file_type'] = 'L1R'
        gatts['saa'] = saa
        gatts['sza'] = sza
        gatts['vaa'] = vaa
        gatts['vza'] = vza
        gatts['raa'] = raa

        ## output name
        oname = '{}_{}'.format(gatts['sensor'], dt.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## new gem handling
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gatts[k] for k in gatts}
        gemo.nc_projection = nc_projection

        ## write position if data are projected
        if dct is not None:
            lon, lat = ac.shared.projection_geo(dct, add_half_pixel = False)
            gemo.write('lon', lon)
            del lon
            gemo.write('lat', lat)
            del lat

        ## printout warning for HJ2B_CCD4
        if sensor == 'HJ2B_CCD4':
            print('Warning: TOA calibration coefficients for HJ2B_CCD4 copied from HJ2B_CCD1.')
            print('Warning: Use results with caution.')

        ## run through bands
        band_indices = [3, 2, 1, 4, 5] ## not correct, tiff labels are incorrect
        band_indices = [1, 2, 3, 4, 5]
        for bi, band in enumerate(bands):
            ## read DN from tiff file, warp if needed
            d = ac.shared.read_band(imagefile, idx = band_indices[bi], sub = sub, rpc_dem = rpc_dem, warp_to = warp_to)
            ## convert to radiance
            d = d.astype(np.float32) * radiometric_calibration[sensor]['scale'][bi] + radiometric_calibration[sensor]['offset'][bi]

            ## output dataset attributes
            ds_att = {'wavelength':rsrd['wave_nm'][band], 'f0': f0d[band]}

            ## apply gains
            if (setu['gains']):
                if len(setu['gains_toa']) == len(bands):
                    cg = float(setu['gains_toa'][bi])
                    ds_att['gain'] = cg
                    if setu['verbosity'] > 2: print('Applying gain {:.5f} for {}'.format(cg, ds))
                    d *= cg

            ## output radiance
            if setu['output_lt']:
                ds = 'Lt_{}'.format(rsrd['wave_name'][band])
                gemo.write(ds, d, ds_att = a)

            ## convert to reflectance
            scale = (np.pi * se_distance ** 2) / (f0d[band] * mus)
            d *= scale

            ## mask negative TOA
            d[d <= 0] = np.nan

            ## write dataset
            ds = 'rhot_{}'.format(rsrd['wave_name'][band])
            gemo.write(ds, d, ds_att = ds_att)
            if setu['verbosity'] > 1: print('Converting bands: Wrote {} ({})'.format(ds, d.shape))
            d = None

        ## write centre position if not determined yet (i.e. unprojected data)
        if (metadata['ProductLevel'] == 'LEVEL1A') & (dct is None):
            gemo.write('lon', centre_lon)
            gemo.write('lat', centre_lat)

        ## close output file
        gemo.close()
        gemo = None

        if ofile is not None:
            if ofile not in ofiles: ofiles.append(ofile)
    return(ofiles, setu)
