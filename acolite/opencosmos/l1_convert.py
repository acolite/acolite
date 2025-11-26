## def l1_convert
## converts OpenCosmos data to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2025-03-20
## modifications: 2025-04-14 (QV) added nc_projection
##                2025-05-20 (QV) test timestamp length
##                2025-05-24 (QV) added scene centre geometry
##                2025-07-22 (QV) added support for ECEF frames

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
        sensor, paths = ac.opencosmos.bundle_test(bundle)

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

        session = list(metadata['Sessions'].keys())[0]
        platform_time = metadata['Sessions'][session]['TimeSync'][0]['PlatformTime']

        if sensor in ['OpenCosmos_Hammer']:
            prj_band = 'HS0'
            bands_fwhm = metadata['Sessions'][session]['ImagerConfiguration']['BandSetup']
            bands_centre = metadata['Sessions'][session]['ImagerConfiguration']['BandCWL']
            nbands = metadata['Sessions'][session]['ImagerConfiguration']['SpectralBands']
            bands_list = list(paths['bands'].keys())
            ## create RSR
            rsr = ac.shared.rsr_hyper(bands_centre, bands_fwhm, names = bands_list, step = 0.1)
            rsrd = ac.shared.rsr_dict(rsrd = {sensor:{'rsr':rsr}})
        else:
            prj_band = 'NIR'
            ## read RSR
            rsrd = ac.shared.rsr_dict(sensor)
        del metadata

        ## test band files and make bands dataset
        bands = {}
        for b in rsrd[sensor]['rsr_bands']:
            if b in paths['bands']:
                if os.path.isfile(paths['bands'][b]):
                    #bands[b] = {'path': paths['bands'][b]}
                    cwave = rsrd[sensor]['wave_nm'][b]
                    swave = '{:.0f}'.format(cwave)
                    bands[b]= {'wave':cwave, 'wavelength':cwave, 'wave_mu':cwave/1000.,
                                   'wave_name':'{:.0f}'.format(cwave), }
                    #           'rsr': rsrd[gatts['sensor']]['rsr'][b],'f0': f0d[b]}
                    bands[b]['path'] = paths['bands'][b]

            if b not in bands:
                print('Band {} not found at {}'.format(b, bundle))

        ## parse date
        time_scale = len('{}'.format(platform_time)) % 10
        time = datetime.datetime.fromtimestamp(platform_time/ 10**time_scale, tz=datetime.timezone.utc)
        doy = time.strftime('%j')
        se_distance = ac.shared.distance_se(doy)

        ## set up projection
        warp_to, dct_prj, sub = None, None, None
        try:
            ## get projection from image
            dct = ac.shared.projection_read(bands[prj_band]['path'])
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
        print('Computing sun geometry')
        spos = ac.shared.sun_position(time, lon, lat)
        sza = spos['zenith'] * 1.0
        saa = spos['azimuth'] * 1.0
        del spos

        ## read ancillary data for satellite position
        print('Computing view geometry')
        with open(paths['ancillary'], 'r', encoding = 'utf-8') as f:
            ancillary = json.load(f)

        ## earth centred earth fixed if GNSS lock achieved
        ecef = ancillary['extrinsics']['gnss_lock']

        ## compute viewing angle at mid point
        mid1 = int(lon.shape[0]/2)
        mid2 = int(lon.shape[1]/2)
        obs_lon = lon[mid1, mid2]
        obs_lat = lat[mid1, mid2]
        obs_alt = 0.

        ## get altitude from DEM
        dem = None
        if setu['dem_altitude']:
            dem = ac.dem.dem_lonlat(lon, lat, source = setu['dem_source'])
            if dem is not None: obs_alt = dem[mid1, mid2]

        ## satellite positions
        satellite_positions = []
        for hi, h in enumerate(ancillary['extrinsics']['hist']):
            if 'position' in h:
                dt =  datetime.datetime.fromtimestamp(h['time_s'], tz=datetime.timezone.utc)
                if h['position']['frame'] == 2: ## frame 2 is ECEF
                    satpos = ac.shared.position.ecef_geometry(
                                        h['position']['x_m']/1000,
                                        h['position']['y_m']/1000,
                                        h['position']['z_m']/1000,
                                        obs_lon, obs_lat, obs_alt)
                elif h['position']['frame'] == 3: ## frame 3 is ECI
                    satpos = ac.shared.position.eci_geometry(
                                        h['position']['x_m']/1000,
                                        h['position']['y_m']/1000,
                                        h['position']['z_m']/1000,
                                        obs_lon, obs_lat, obs_alt, dt)
                else:
                    print('Frame not recognised: {}'.format(h))
                    continue

                ## append position
                satellite_positions.append([h['time_s'], dt, satpos])

        ## find closest element in sat time
        sat_time = np.asarray([s[0] for s in satellite_positions])
        idx = np.argsort(np.abs((sat_time - (platform_time))))[0]
        vza = satellite_positions[idx][2]['zenith']
        vaa = satellite_positions[idx][2]['azimuth']

        raa = np.abs(saa - vaa)
        raa[raa>180]= np.abs(360-raa[raa>180])

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

        ## add band config for HAMMER
        if sensor in ['OpenCosmos_Hammer']:
            gatts['band_widths'] = bands_fwhm
            gatts['band_waves'] = bands_centre

        ## output file name
        oname = '{}_{}'.format(gatts['sensor'], time.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## DN to reflectance
        reflectance_scale_factor = 0.0001 ## email Melina 2025-03-14

        ## set up output gem
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gatts[k] for k in gatts}
        gemo.nc_projection = nc_projection

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

        ## output geolocation
        if setu['output_geometry']:
            if setu['verbosity'] > 1: print('Writing geometry')
            gemo.write('sza', sza)
            gemo.write('saa', saa)
            gemo.write('vza', vza)
            gemo.write('vaa', vaa)
            gemo.write('raa', raa)
            del sza, saa, vza, vaa, raa

        ## run through bands
        for bi, b in enumerate(bands):
            if setu['verbosity'] > 2: print('Computing rhot_{} for {}'.format(bands[b]['wave_name'], gatts['oname']))
            ds_att = {k: bands[b][k] for k in bands[b] if k not in ['rsr']}
            if 'path' in ds_att: del ds_att['path']

            wave = bands[b]['wavelength']
            ds = 'rhot_{:.0f}'.format(wave)

            ## read data and compute mask
            cdata = ac.shared.read_band(bands[b]['path'], sub = sub, warp_to = warp_to)
            mask = cdata == 0

            cdata = cdata.astype(np.float32) * reflectance_scale_factor
            cdata[mask] = np.nan

            gemo.write('rhot_{}'.format(bands[b]['wave_name']), cdata, ds_att = ds_att)
            cdata = None
        gemo.close()

        ofiles.append(ofile)
    return(ofiles, setu)
