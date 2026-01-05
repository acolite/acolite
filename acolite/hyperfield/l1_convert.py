## def l1_convert
## converts Hyperfield data to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2026-01-06
## modifications:

def l1_convert(inputfile, output = None, settings = None):
    import os, json
    import numpy as np
    import dateutil.parser, datetime
    import acolite as ac

    from osgeo import gdal
    gdal.UseExceptions()

    ## get run/user/sensor settings
    setu = ac.acolite.settings.merge(sensor = None, settings = settings)

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if setu['verbosity'] > 1: print('Starting conversion of {} scene{}'.format(nscenes, '' if nscenes==1 else 's'))

    ## list to store output files
    ofiles = []

    ## run through inputfiles
    for bi, bundle in enumerate(inputfile):
        ## identify sensor and file paths
        sensor_, paths = ac.hyperfield.bundle_test(bundle)
        aliases = {'hyperfield1a': 'Hyperfield-1A',
                   'hyperfield1b': 'Hyperfield-1B',}
        if sensor_ in aliases:
            sensor = aliases[sensor_]
        else:
            print('Sensor {} not configured'.format(sensor_))
            continue

        ## update settings
        setu = ac.acolite.settings.merge(sensor = sensor, settings = settings)

        ## set output directory
        if output is None: output = setu['output']
        if output is None: output = os.path.dirname(file)

        ## read F0
        #if setu['verbosity'] > 2: print('Using F0 from solar_irradiance_reference={}'.format(setu['solar_irradiance_reference']))
        #f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        #f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], rsrd[gatts['sensor']]['rsr'])

        ## read metadata to get timestamp (and HAMMER bands config)
        with open(paths['metadata'], 'r', encoding = 'utf-8') as f:
            metadata = json.load(f)

        if metadata['satellite']['name'] != sensor.lower():
            print('Sensor {} does not match metadata {} not configured'.format(sensor, metadata['satellite']))
            continue

        if metadata['image']['measured_quantity_name'] != 'TOA_REFLECTANCE':
            print('Measured quantity {} not supported'.format(metadata['image']['measured_quantity_name']))
            continue

        ## image info
        isodate = metadata['image']['acquired_on']
        time = dateutil.parser.parse(isodate)
        doy = time.strftime('%j')
        se_distance = ac.shared.distance_se(doy)

        ## extract bands info
        bands_fwhm = []
        bands_centre = []
        nbands = len(metadata['image']['bands'])
        for b in metadata['image']['bands']:
            bands_fwhm.append(b['bandwidth'][0])
            bands_centre.append(b['wavelength'][0])
        ## create rsr
        rsr = ac.shared.rsr_hyper(bands_centre, bands_fwhm, step = 0.1)
        rsrd = ac.shared.rsr_dict(rsrd = {sensor:{'rsr':rsr}})
        ## create bands dict
        bands = {}
        for b in rsrd[sensor]['rsr_bands']:
            cwave = rsrd[sensor]['wave_nm'][b]
            swave = '{:.0f}'.format(cwave)
            bands[b]= {'wave':cwave, 'wavelength':cwave, 'wave_mu':cwave/1000.,
                                     'wave_name':'{:.0f}'.format(cwave), }
            for k in metadata['image']['bands'][b]:
                if k not in bands[b]:
                    bands[b][k] = metadata['image']['bands'][b][k]

        ## set up projection
        warp_to, dct_prj, sub = None, None, None
        try:
            ## get projection from image
            dct = ac.shared.projection_read(paths['L1C'])
        except:
            if setu['verbosity'] > 1: print('Could not determine image projection')
            dct = None

        ## find crop
        if (setu['limit'] is not None) and (dct is not None):
            dct_sub = ac.shared.projection_sub(dct, setu['limit'])
            if dct_sub['out_lon']:
                if setu['verbosity'] > 1: print('Longitude limits outside {}'.format(bundle))
                continue
            if dct_sub['out_lat']:
                if setu['verbosity'] > 1: print('Latitude limits outside {}'.format(bundle))
                continue
            sub = dct_sub['sub']

        nc_projection = None
        if dct is not None:
            if sub is None:
                dct_prj = {k:dct[k] for k in dct}
            else:
                dct_prj = {k:dct_sub[k] for k in dct_sub}
                xyr = [min(dct_prj['xrange']),
                       min(dct_prj['yrange']),
                       max(dct_prj['xrange']),
                       max(dct_prj['yrange']),
                       dct_prj['proj4_string']]
                res_method = 'near'
                warp_to = (dct_prj['proj4_string'], xyr, dct_prj['pixel_size'][0],dct_prj['pixel_size'][1], res_method)
            if setu['netcdf_projection']:
                nc_projection = ac.shared.projection_netcdf(dct_prj, add_half_pixel=False)

        ## compute lat/lon
        lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel = False)

        ## compute sun angles
        # print('Computing sun geometry')
        # spos = ac.shared.sun_position(time, lon, lat)
        # sza = spos['zenith'] * 1.0
        # saa = spos['azimuth'] * 1.0
        # del spos

        ## scene centre and average geometry
        sza = metadata['image']['local_solar_zenith_angle'][0]
        saa = metadata['image']['local_solar_azimuth_angle'][0]
        vza_list = [bands[b]['viewing_zenith_angle'][0] for b in bands]
        vaa_list = [bands[b]['viewing_azimuth_angle'][0] for b in bands]
        vza = np.nanmean(vza_list)
        vaa = np.nanmean(vaa_list)
        raa = np.abs(saa - vaa)
        #raa[raa>180]= np.abs(360-raa[raa>180])
        if raa > 180: raa = np.abs(360-raa)

        #del metadata

        ## output gatts
        gatts = {}
        gatts['sensor'] = sensor
        gatts['isodate'] = time.isoformat()
        gatts['acolite_file_type'] = 'L1R'
        gatts['doy'] = doy
        gatts['se_distance'] = se_distance
        gatts['acolite_file_type'] =  'L1R'
        gatts['sza'] = np.nanmean(sza)
        gatts['saa'] = np.nanmean(saa)
        gatts['vza'] = np.nanmean(vza)
        gatts['vaa'] = np.nanmean(vaa)
        gatts['raa'] = np.nanmean(raa)

        ## add band config
        gatts['band_widths'] = bands_fwhm
        gatts['band_waves'] = bands_centre

        ## output file name
        oname = '{}_{}'.format(gatts['sensor'], time.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## set up output gem
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gatts[k] for k in gatts}
        gemo.nc_projection = nc_projection

        dem = None

        ## output geolocation
        if setu['output_geolocation']:
            if setu['verbosity'] > 1: print('Writing lat/lon')
            gemo.write('lat', lat)
            lat = None
            gemo.write('lon', lon)
            lon = None
            if dem is not None:
                gemo.write('dem', dem)
                dem = None

        # ## output geolocation
        # if setu['output_geometry']:
        #     if setu['verbosity'] > 1: print('Writing geometry')
        #     gemo.write('sza', sza)
        #     gemo.write('saa', saa)
        #     gemo.write('vza', vza)
        #     gemo.write('vaa', vaa)
        #     gemo.write('raa', raa)
        #     del sza, saa, vza, vaa, raa


        ## run through bands
        for bi, b in enumerate(bands):
            if setu['verbosity'] > 2: print('Computing rhot_{} for {}'.format(bands[b]['wave_name'], gatts['oname']))
            ds_att = {k: bands[b][k] for k in bands[b] if k not in ['rsr']}
            if 'path' in ds_att: del ds_att['path']

            wave = bands[b]['wavelength']
            ds = 'rhot_{:.0f}'.format(wave)

            ## read data and compute mask
            cdata = ac.shared.read_band(paths['L1C'], idx = bands[b]['index'] + 1, sub = sub, warp_to = warp_to)
            mask = cdata == 0
            cdata = cdata.astype(np.float32)
            if bands[b]['offset'] != 0.0: cdata += bands[b]['offset']
            if bands[b]['scale'] != 1.0: cdata *= bands[b]['scale']

            #cdata /= bands[b]['toa_radiance_to_reflectance_factor']
            cdata[mask] = np.nan

            gemo.write('rhot_{}'.format(bands[b]['wave_name']), cdata, ds_att = ds_att)
            cdata = None
        gemo.close()

        ofiles.append(ofile)
    return(ofiles, setu)
