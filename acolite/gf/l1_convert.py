## def l1_convert
## converts GF6 bundle to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-08-09
## modifications: 2021-11-20 (QV) reproject file if projection not recognised
##                2021-12-31 (QV) new handling of settings
##                2022-01-04 (QV) added netcdf compression
##                2022-11-22 (QV) added GF1 WFV1-4
##                2022-12-10 (QV) changed bias to 0.0
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-17 (QV) use new gem NetCDF handling
##                2025-01-28 (QV) switch to LinearNDInterpolator, added meshgrid
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming

def l1_convert(inputfile, output = None, settings = None):
    import numpy as np
    import scipy.interpolate

    import datetime, dateutil.parser, os
    import acolite as ac
    from osgeo import gdal

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

    ## from Mona Allam PDF
    ## 240A9946ECED64C3CB4D56B8A8764F35
    dn_scaling = {}
    dn_scaling['GF1'] = {}
    dn_scaling['GF1']['WFV1'] = {'Blue':0.19319, 'Green': 0.16041, 'Red': 0.12796, 'NIR': 0.13405}
    dn_scaling['GF1']['WFV2'] = {'Blue':0.2057, 'Green': 0.1648, 'Red': 0.1260, 'NIR': 0.1187}
    dn_scaling['GF1']['WFV3'] = {'Blue':0.2106, 'Green': 0.1825, 'Red': 0.1346, 'NIR': 0.1187}
    dn_scaling['GF1']['WFV4'] = {'Blue':0.2522, 'Green': 0.2029, 'Red': 0.1528, 'NIR': 0.1031}

    dn_scaling['GF1B'] = {'PMS': {'PAN': 0.0687, 'MS1': 0.0757, 'MS2': 0.0618, 'MS3': 0.0545, 'MS4': 0.0572}}
    dn_scaling['GF1C'] = {'PMS': {'PAN': 0.0709, 'MS1': 0.0758, 'MS2': 0.0657, 'MS3': 0.0543, 'MS4': 0.0564}}
    dn_scaling['GF1D'] = {'PMS': {'PAN': 0.0715, 'MS1': 0.0738, 'MS2': 0.0656, 'MS3': 0.0590, 'MS4': 0.0585}}

    dn_scaling['GF6'] = {'WFV': {'B1': 0.0675, 'B2': 0.0552, 'B3': 0.0513, 'B4': 0.0314,
                                 'B5': 0.0519, 'B6': 0.0454, 'B7': 0.0718, 'B8': 0.0596},
                          'PMS': {'PAN': 0.0537, 'MS1': 0.082, 'MS2': 0.0645, 'MS3': 0.0489, 'MS4': 0.0286}}

    ## bias should be 0 and not 0.2
    ## https://github.com/acolite/acolite/issues/53
    dn_bias = 0.0

    ofiles = []
    for bundle in inputfile:
        tiles, metafile = ac.gf.bundle_test(bundle)
        meta = ac.gf.metadata(metafile)
        if meta['SatelliteID'] not in  ['GF1', 'GF1D', 'GF6']: continue
        sensor = '{}_{}'.format(meta['SatelliteID'], meta['SensorID'])

        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults

        verbosity = setu['verbosity']
        if output is None: output = setu['output']
        if output is None: output = os.path.dirname(bundle)

        ## get F0 for radiance -> reflectance computation
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])

        print('Processing {}'.format(bundle))

        ## parse data
        dtime = dateutil.parser.parse(meta['CenterTime'])
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()

        ## figure out suitable UTM zone
        if 'CenterLongitude' in meta:
            clon = float(meta['CenterLongitude'])
        else:
            clon = 0.5*(float(meta['BottomLeftLongitude'])+float(meta['TopRightLongitude']))
        if 'CenterLatitude' in meta:
            clat = float(meta['CenterLatitude'])
        else:
            clat = 0.5*(float(meta['BottomLeftLatitude'])+float(meta['TopRightLatitude']))

        utm_zone = (int(1+(clon+180.0)/6.0))
        north = clat > 0.0
        epsg = 'EPSG:32{}{}'.format('6' if north else '7', utm_zone)

        ## output attributes
        gatts = {}
        gatts['sza'] = float(meta['SolarZenith'])
        gatts['vza'] = float(meta['SatelliteZenith'])
        gatts['saa'] = float(meta['SolarAzimuth'])
        gatts['vaa'] = float(meta['SatelliteAzimuth'])

        if 'raa' not in gatts:
            raa_ave = abs(gatts['saa'] - gatts['vaa'])
            while raa_ave >= 180: raa_ave = abs(raa_ave-360)
            gatts['raa'] = raa_ave

        gatts['satellite'] = meta['SatelliteID']
        gatts['sensor'] = sensor
        gatts['isodate'] = isodate
        gatts['se_distance'] = se_distance
        gatts['doy'] = doy
        gatts['acolite_file_type'] = 'L1R'

        ##
        sensor = gatts['sensor']
        rsrd = ac.shared.rsr_dict(sensor)[sensor]
        band_names = rsrd['rsr_bands']
        f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], rsrd['rsr'])

        mu0 = np.cos(gatts['sza']*(np.pi/180))
        muv = np.cos(gatts['vza']*(np.pi/180))

        bands = {}
        for bi, b in enumerate(band_names):
            bands[b] = {}
            for k in ['wave_mu', 'wave_nm', 'wave_name']:
                bands[b][k] = rsrd[k][b]
            bands[b]['f0'] = f0d[b]

            if sensor in ['GF1_WFV1', 'GF1_WFV2', 'GF1_WFV3', 'GF1_WFV4', 'GF1D_PMS', 'GF6_PMS']:
                bands[b]['index'] = int(bi)+1
            if sensor in ['GF6_WFV']:
                bands[b]['index'] = int(b[1])

        ## order bands
        if sensor in ['GF1D_PMS', 'GF6_PMS']:
            idx = np.argsort([bands[b]['wave_name'] for b in bands if b not in ['PAN']])
        if sensor in ['GF1_WFV1', 'GF1_WFV2', 'GF1_WFV3', 'GF1_WFV4', 'GF6_WFV']:
            idx = np.argsort([bands[b]['wave_name'] for b in bands])
        bands_sorted = [band_names[i] for i in idx]

        ## image crop
        if setu['limit'] is None: sub = None

        ## track tile offsets
        x_off = 0
        y_off = 0

        ## run through tiles
        if verbosity > 1: print('Running through {} {} image tiles'.format(len(tiles), sensor))
        for ti, image_file in enumerate(tiles):
            bn = os.path.basename(image_file)
            try:
                ctile = os.path.splitext(bn)[0].split('-')[1]
            except:
                ctile = meta['ProductID']

            #if ('PMS' in bn) & ('PAN' in bn): continue
            if ctile.upper() == 'PAN': continue
            oname  = '{}_{}_{}'.format(gatts['sensor'],  dtime.strftime('%Y_%m_%d_%H_%M_%S'), ctile)
            if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
            ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
            gatts['oname'] = oname
            gatts['ofile'] = ofile

            ## output file - one per tile
            gemo = ac.gem.gem(ofile, new = True)
            gemo.gatts = {k: gatts[k] for k in gatts}
            #gemo.nc_projection = nc_projection

            if verbosity > 1: print('Running {} tile {}/{}'.format(sensor, ti+1, len(tiles)))

            ## identify projection
            try:
                prj = ac.shared.projection_read(image_file)
                if verbosity > 1: print('Could determine projection from {}'.format(image_file))
            except:
                prj = None

            ## reproject file to UTM if projection not read succesfully
            rpr_file = None
            if setu['gf_reproject_to_utm']:
                if prj is None:
                    rpr_file = '{}/{}'.format(ac.config['scratch_dir'], bn.replace('.tiff', '_reprojected.tif'))
                    if not os.path.exists(rpr_file):
                        if verbosity > 1: print('Reprojecting {} to {}'.format(image_file, epsg))
                        if verbosity > 1: print('Target file {}'.format(rpr_file))
                        ## scratch directory
                        if not os.path.exists(ac.config['scratch_dir']): os.makedirs(ac.config['scratch_dir'])
                        ## reproject and close dataset
                        ds = gdal.Warp(rpr_file, image_file, dstSRS=epsg)
                        ds = None
                    ## get prj from new file
                    if os.path.exists(rpr_file): prj = ac.shared.projection_read(rpr_file)

            ## add projection keys to gatts
            if prj is not None:
                pkeys = ['xrange', 'yrange', 'proj4_string', 'pixel_size', 'zone']
                for k in pkeys:
                    if k in prj: gatts[k] = prj[k]

                ## compute geolocation
                if True:
                    if verbosity > 1: print('Computing latitude/longitude')
                    lon, lat = ac.shared.projection_geo(prj, add_half_pixel=True)
                    gemo.write('lon', lon)
                    lon = None
                    if verbosity > 1: print('Wrote lon')
                    gemo.write('lat', lat)
                    lat = None
                    if verbosity > 1: print('Wrote lat')

            ## run through bands
            for bi, b in enumerate(bands_sorted):
                if b == 'PAN': continue
                print('Computing rhot_{} for {}'.format(bands[b]['wave_name'], gatts['oname']))
                print(b, bi, dn_scaling[gatts['satellite']][gatts['sensor'].split('_')[1]][b])

                ## read data
                if rpr_file is None:
                    cdata_radiance = ac.shared.read_band(image_file, bands[b]['index'], sub=sub)
                else:
                    cdata_radiance = ac.shared.read_band(rpr_file, bands[b]['index'], sub=sub)
                data_shape = cdata_radiance.shape

                ## compute radiance
                cdata_radiance = cdata_radiance.astype(np.float32) * dn_scaling[gatts['satellite']][gatts['sensor'].split('_')[1]][b]
                cdata_radiance += dn_bias

                if setu['output_lt']:
                    ## write toa radiance
                    gemo.write('Lt_{}'.format(bands[b]['wave_name']), cdata_radiance, ds_att = bands[b])

                ## compute reflectance
                cdata = cdata_radiance * (np.pi * gatts['se_distance'] * gatts['se_distance']) / (bands[b]['f0'] * mu0)
                cdata_radiance = None
                gemo.write('rhot_{}'.format(bands[b]['wave_name']), cdata, ds_att = bands[b])
                cdata = None

            ## old geolocation
            if (rpr_file is None) and (prj is None):
                nx, ny = int(meta['WidthInPixels']), int(meta['HeightInPixels'])
                tllat, tllon = float(meta['TopLeftLatitude']), float(meta['TopLeftLongitude'])
                trlat, trlon = float(meta['TopRightLatitude']), float(meta['TopRightLongitude'])
                bllat, bllon = float(meta['BottomLeftLatitude']), float(meta['BottomLeftLongitude'])
                brlat, brlon = float(meta['BottomRightLatitude']), float(meta['BottomRightLongitude'])
                #clat, clon = float(meta['CenterLatitude']), float(meta['CenterLongitude'])

                ## get vertex image location
                pcol = [0, nx, nx, 0]
                prow = [0, 0, ny, ny]

                ## get vertex geolocation
                plon = [tllon, trlon, brlon, bllon]
                plat = [tllat, trlat, brlat, bllat]

                ## set up interpolator
                zlon = scipy.interpolate.LinearNDInterpolator((pcol, prow), plon)
                zlat = scipy.interpolate.LinearNDInterpolator((pcol, prow), plat)

                ## pixel coordinate limits
                if sensor in ['GF1_WFV1', 'GF1_WFV2', 'GF1_WFV3', 'GF1_WFV4','GF1D_PMS', 'GF6_PMS']:
                    x0, y0 = 0, 0
                    ns, nl = nx, ny
                if sensor == 'GF6_WFV':
                    x0 = x_off
                    y0 = y_off
                    ns, nl = data_shape[1], data_shape[0]
                    x_off += data_shape[1]

                x = np.arange(x0, x0+ns, 1)
                y = np.arange(y0, y0+nl, 1)
                X, Y = np.meshgrid(x, y)

                print('Computing lon')
                gemo.write('lon', zlon(X, Y))

                print('Computing lat')
                gemo.write('lat', zlat(X, Y))
                del X, Y

            ## remove reprojected file
            if (setu['clear_scratch']) & (rpr_file is not None):
                if os.path.exists(rpr_file):
                    os.remove(rpr_file)

            ## close output file
            gemo.close()

            ## add current tile to outputs
            ofiles.append(ofile)

    ## remove scratch directory
    if os.path.exists(ac.config['scratch_dir']):
        if (setu['clear_scratch']) & (len(os.listdir(ac.config['scratch_dir'])) == 0):
            os.rmdir(ac.config['scratch_dir'])

    return(ofiles, setu)
