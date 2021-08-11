## def l1_convert
## converts GF6 bundle to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-08-09

def l1_convert(inputfile, output = None, limit = None, verbosity=0, vname = '', output_lt=False):
    import numpy as np
    from scipy.interpolate import interp2d

    import datetime, dateutil.parser, os
    import acolite as ac

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
    dn_scaling = {'WFV': {'B1': 0.0675, 'B2': 0.0552, 'B3': 0.0513, 'B4': 0.0314,
                          'B5': 0.0519, 'B6': 0.0454, 'B7': 0.0718, 'B8': 0.0596},
                  'PMS': {'PAN': 0.0537, 'MS1': 0.082, 'MS2': 0.0645, 'MS3': 0.0489, 'MS4': 0.0286}}
    dn_bias = 0.2

    ofiles = []
    for bundle in inputfile:
        tiles, metafile = ac.gf6.bundle_test(bundle)
        meta = ac.gf6.metadata(metafile)
        if meta['SatelliteID'] != 'GF6': continue

        print('Processing {}'.format(bundle))

        ## parse data
        dtime = dateutil.parser.parse(meta['CenterTime'])
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()

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
            if sensor == 'GF6_WFV':
                bands[b]['index'] = int(b[1])
            if sensor == 'GF6_PMS':
                bands[b]['index'] = int(bi)+1

        ## order bands
        if sensor == 'GF6_WFV':
            idx = np.argsort([bands[b]['wave_name'] for b in bands])
        if sensor == 'GF6_PMS':
            idx = np.argsort([bands[b]['wave_name'] for b in bands if b not in ['PAN']])
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

            ## run through bands
            for bi, b in enumerate(bands_sorted):
                if b == 'PAN': continue
                print('Computing rhot_{} for {}'.format(bands[b]['wave_name'], gatts['obase']))
                print(b, bi, dn_scaling[gatts['sensor'].split('_')[1]][b])

                ## read data
                cdata_radiance = ac.shared.read_band(image_file, bands[b]['index'], sub=sub)
                data_shape = cdata_radiance.shape

                ## compute radiance
                cdata_radiance = cdata_radiance.astype(np.float32) * dn_scaling[gatts['sensor'].split('_')[1]][b]
                cdata_radiance += dn_bias

                if output_lt:
                    ## write toa radiance
                    ac.output.nc_write(ofile, 'Lt_{}'.format(bands[b]['wave_name']), cdata_radiance,
                                        attributes = gatts, dataset_attributes = bands[b], new = new)
                    new = False

                ## compute reflectance
                cdata = cdata_radiance * (np.pi * gatts['se_distance'] * gatts['se_distance']) / (bands[b]['f0'] * mu0)
                cdata_radiance = None

                ac.output.nc_write(ofile, 'rhot_{}'.format(bands[b]['wave_name']), cdata,\
                                        attributes = gatts, dataset_attributes = bands[b], new = new)
                cdata = None
                new = False

            ## geolocation
            if True:
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
                if sensor == 'GF6_PMS':
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
                ac.output.nc_write(ofile, 'lat', lat, attributes = gatts, new = new)
                lat = None

                print('Computing lon')
                lon = zlon(x, y)
                ac.output.nc_write(ofile, 'lon', lon)
                lon = None
                new = False


            ## add current tile to outputs
            ofiles.append(ofile)
    return(ofiles)
