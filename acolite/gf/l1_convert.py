## def l1_convert
## converts GF6 bundle to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-08-09
## modifications: 2021-11-20 (QV) reproject file if projection not recognised
##                2021-12-31 (QV) new handling of settings
##                2022-01-04 (QV) added netcdf compression

def l1_convert(inputfile, output = None, settings = {}, verbosity=5):
    import numpy as np
    from scipy.interpolate import interp2d

    import datetime, dateutil.parser, os
    import acolite as ac
    from osgeo import gdal
    #import subprocess

    if 'verbosity' in settings: verbosity = settings['verbosity']

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ## get F0 for radiance -> reflectance computation
    f0 = ac.shared.f0_get()

    ## from Mona Allam PDF
    ## 240A9946ECED64C3CB4D56B8A8764F35
    dn_scaling = {}
    dn_scaling['GF1B'] = {'PMS': {'PAN': 0.0687, 'MS1': 0.0757, 'MS2': 0.0618, 'MS3': 0.0545, 'MS4': 0.0572}}
    dn_scaling['GF1C'] = {'PMS': {'PAN': 0.0709, 'MS1': 0.0758, 'MS2': 0.0657, 'MS3': 0.0543, 'MS4': 0.0564}}
    dn_scaling['GF1D'] = {'PMS': {'PAN': 0.0715, 'MS1': 0.0738, 'MS2': 0.0656, 'MS3': 0.0590, 'MS4': 0.0585}}

    dn_scaling['GF6'] = {'WFV': {'B1': 0.0675, 'B2': 0.0552, 'B3': 0.0513, 'B4': 0.0314,
                                 'B5': 0.0519, 'B6': 0.0454, 'B7': 0.0718, 'B8': 0.0596},
                          'PMS': {'PAN': 0.0537, 'MS1': 0.082, 'MS2': 0.0645, 'MS3': 0.0489, 'MS4': 0.0286}}
    dn_bias = 0.2

    ofiles = []
    for bundle in inputfile:
        tiles, metafile = ac.gf.bundle_test(bundle)
        meta = ac.gf.metadata(metafile)
        if meta['SatelliteID'] not in  ['GF1D', 'GF6']: continue

        ## sensor settings
        setu = ac.acolite.settings.parse(meta['SatelliteID'], settings=settings)
        verbosity = setu['verbosity']
        ## get other settings
        limit = setu['limit']
        output_lt = setu['output_lt']
        reproject_to_utm = setu['gf_reproject_to_utm']
        clear_scratch = setu['clear_scratch']
        vname = setu['region_name']
        gains = setu['gains']
        gains_toa = setu['gains_toa']
        if output is None: output = setu['output']

        print('Processing {}'.format(bundle))

        ## parse data
        dtime = dateutil.parser.parse(meta['CenterTime'])
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()

        ## figure out suitable UTM zone
        utm_zone = (int(1+(float(meta['CenterLongitude'])+180.0)/6.0))
        north = float(meta['CenterLatitude']) > 0.0
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
        gatts['sensor'] = '{}_{}'.format(meta['SatelliteID'], meta['SensorID'])
        gatts['isodate'] = isodate
        gatts['se_distance'] = se_distance
        gatts['doy'] = doy

        obase  = '{}_{}_L1R'.format(gatts['sensor'],  dtime.strftime('%Y_%m_%d_%H_%M_%S'))

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

            if sensor in ['GF1D_PMS', 'GF6_PMS']:
                bands[b]['index'] = int(bi)+1
            if sensor == 'GF6_WFV':
                bands[b]['index'] = int(b[1])

        ## order bands
        if sensor in ['GF1D_PMS', 'GF6_PMS']:
            idx = np.argsort([bands[b]['wave_name'] for b in bands if b not in ['PAN']])
        if sensor == 'GF6_WFV':
            idx = np.argsort([bands[b]['wave_name'] for b in bands])
        bands_sorted = [band_names[i] for i in idx]

        ## image crop
        if limit is None: sub = None

        ## output file
        if output is None:
            odir = os.path.dirname(bundle)
        else:
            odir = '{}'.format(output)
        if not os.path.exists(odir): os.makedirs(odir)

        ## track tile offsets
        x_off = 0
        y_off = 0

        ## run through tiles
        for ti, image_file in enumerate(tiles):
            new = True
            bn = os.path.basename(image_file)
            ctile = os.path.splitext(bn)[0].split('-')[1]
            #if ('PMS' in bn) & ('PAN' in bn): continue
            if ctile.upper() == 'PAN': continue
            gatts['obase'] = obase + '_{}'.format(ctile)
            ofile = '{}/{}.nc'.format(odir, gatts['obase'])

            ## identify projection
            try:
                prj = ac.shared.projection_read(image_file)
            except:
                prj = None

            ## reproject file to UTM if projection not read succesfully
            rpr_file = None
            if reproject_to_utm:
                if prj is None:
                    rpr_file = '{}/{}'.format(ac.config['scratch_dir'], bn.replace('.tiff', '_reprojected.tif'))
                    if not os.path.exists(rpr_file):
                        if verbosity > 1: print('Reprojecting {} to {}'.format(image_file, epsg))
                        if verbosity > 1: print('Target file {}'.format(rpr_file))
                        #sp = subprocess.run(' '.join(["gdalwarp", "-overwrite", " -t_srs {} -r cubic ".format(epsg),
                        #                              image_file.replace(' ', '\ '), rpr_file.replace(' ', '\ ')]),
                        #                              shell=True, check=True, stdout=subprocess.PIPE)
                        ## scratch directory
                        if not os.path.exists(ac.config['scratch_dir']):
                            os.makedirs(ac.config['scratch_dir'])
                        ## reproject and close dataset
                        ds = gdal.Warp(rpr_file, image_file, dstSRS=epsg)
                        ds = None
                    ## get prj from new file
                    if os.path.exists(rpr_file):
                        prj = ac.shared.projection_read(rpr_file)

            ## add projection keys to gatts
            if prj is not None:
                pkeys = ['xrange', 'yrange', 'proj4_string', 'pixel_size', 'zone']
                for k in pkeys:
                    if k in prj: gatts[k] = prj[k]

                ## compute geolocation
                if True:
                    if verbosity > 1: print('Computing latitude/longitude')
                    lon, lat = ac.shared.projection_geo(prj, add_half_pixel=True)
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
                    new=False

            ## run through bands
            for bi, b in enumerate(bands_sorted):
                if b == 'PAN': continue
                print('Computing rhot_{} for {}'.format(bands[b]['wave_name'], gatts['obase']))
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

                if output_lt:
                    ## write toa radiance
                    ac.output.nc_write(ofile, 'Lt_{}'.format(bands[b]['wave_name']), cdata_radiance,
                                        attributes = gatts, dataset_attributes = bands[b], new = new,
                                        netcdf_compression=setu['netcdf_compression'],
                                        netcdf_compression_level=setu['netcdf_compression_level'],
                                        netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
                    new = False

                ## compute reflectance
                cdata = cdata_radiance * (np.pi * gatts['se_distance'] * gatts['se_distance']) / (bands[b]['f0'] * mu0)
                cdata_radiance = None

                ac.output.nc_write(ofile, 'rhot_{}'.format(bands[b]['wave_name']), cdata,\
                                        attributes = gatts, dataset_attributes = bands[b], new = new,
                                        netcdf_compression=setu['netcdf_compression'],
                                        netcdf_compression_level=setu['netcdf_compression_level'],
                                        netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
                cdata = None
                new = False

            ## old geolocation
            if rpr_file is None:
                nx, ny = int(meta['WidthInPixels']), int(meta['HeightInPixels'])
                tllat, tllon = float(meta['TopLeftLatitude']), float(meta['TopLeftLongitude'])
                trlat, trlon = float(meta['TopRightLatitude']), float(meta['TopRightLongitude'])
                bllat, bllon = float(meta['BottomLeftLatitude']), float(meta['BottomLeftLongitude'])
                brlat, brlon = float(meta['BottomRightLatitude']), float(meta['BottomRightLongitude'])
                clat, clon = float(meta['CenterLatitude']), float(meta['CenterLongitude'])

                ## get vertex image location
                pcol = [0, nx, nx, 0]
                prow = [0, 0, ny, ny]

                ## get vertex geolocation
                plon = [tllon, trlon, brlon, bllon]
                plat = [tllat, trlat, brlat, bllat]

                ## set up interpolator
                zlon = interp2d(pcol, prow, plon)
                zlat = interp2d(pcol, prow, plat)

                ## pixel coordinate limits
                if sensor in ['GF1D_PMS', 'GF6_PMS']:
                    x0, y0 = 0, 0
                    ns, nl = nx, ny
                if sensor == 'GF6_WFV':
                    x0 = x_off
                    y0 = y_off
                    ns, nl = data_shape[1], data_shape[0]
                    x_off += data_shape[1]

                x = np.arange(x0, x0+ns, 1)
                y = np.arange(y0, y0+nl, 1)

                print('Computing lat')
                lat = zlat(x, y)
                ac.output.nc_write(ofile, 'lat', lat, attributes = gatts, new = new,
                                    netcdf_compression=setu['netcdf_compression'],
                                    netcdf_compression_level=setu['netcdf_compression_level'])
                lat = None

                print('Computing lon')
                lon = zlon(x, y)
                ac.output.nc_write(ofile, 'lon', lon,
                                    netcdf_compression=setu['netcdf_compression'],
                                    netcdf_compression_level=setu['netcdf_compression_level'])
                lon = None
                new = False

            ## remove reprojected file
            if (clear_scratch) & (rpr_file is not None):
                if os.path.exists(rpr_file):
                    os.remove(rpr_file)

            ## add current tile to outputs
            ofiles.append(ofile)

    ## remove scratch directory
    if os.path.exists(ac.config['scratch_dir']):
        if (clear_scratch) & (len(os.listdir(ac.config['scratch_dir'])) == 0):
            os.rmdir(ac.config['scratch_dir'])

    return(ofiles, setu)
