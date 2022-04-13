## def l1_convert
## converts FORMOSAT5 bundle to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2022-04-12
## modifications:

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


    ofiles = []
    for bundle in inputfile:
        tiles, metafile = ac.formosat.bundle_test(bundle)
        meta = ac.formosat.metadata(metafile)
        if meta['sensor'] not in ['FORMOSAT5_RSI']: continue

        ## sensor settings
        setu = ac.acolite.settings.parse(meta['sensor'], settings=settings)
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
        dtime = dateutil.parser.parse(meta['isodate'])
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()

        ## output attributes
        gatts = {}
        gatts['sza'] = float(meta['sza'])
        gatts['vza'] = float(meta['vza'])
        gatts['saa'] = float(meta['saa'])
        gatts['vaa'] = float(meta['vaa'])

        if 'raa' not in gatts:
            raa_ave = abs(gatts['saa'] - gatts['vaa'])
            while raa_ave >= 180: raa_ave = abs(raa_ave-360)
            gatts['raa'] = raa_ave

        gatts['satellite'] = meta['satellite']
        gatts['sensor'] = meta['sensor']
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
        band_indices = {'Blue':3, 'Green':2, 'Red':1, 'NIR':4, 'PAN': -1}

        for bi, b in enumerate(band_names):
            bands[b] = {}
            for k in ['wave_mu', 'wave_nm', 'wave_name']:
                bands[b][k] = rsrd[k][b]
            bands[b]['f0'] = f0d[b]
            bands[b]['index'] = band_indices[b]

        ## image crop
        if limit is None: sub = None

        ## output file
        if output is None:
            odir = os.path.dirname(bundle)
        else:
            odir = '{}'.format(output)
        if not os.path.exists(odir): os.makedirs(odir)

        ## run through tiles
        for ti, image_file in enumerate(tiles):
            new = True
            gatts['obase'] = obase
            ofile = '{}/{}.nc'.format(odir, gatts['obase'])

            ## identify projection
            try:
                prj = ac.shared.projection_read(image_file)
            except:
                prj = None

            ## add projection keys to gatts
            if prj is not None:
                pkeys = ['xrange', 'yrange', 'proj4_string', 'pixel_size', 'zone']
                for k in pkeys:
                    if k in prj: gatts[k] = prj[k]

                ## compute geolocation
                if setu['output_geolocation']:
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
            for bi, b in enumerate(band_names):
                if b == 'PAN': continue
                gain = meta['{}_PHYSICAL_GAIN'.format(b)]
                bias = meta['{}_PHYSICAL_BIAS'.format(b)]
                print('Computing rhot_{} for {}'.format(bands[b]['wave_name'], gatts['obase']))
                print(b, bi, gain, bias)

                ## read data
                cdata_radiance = ac.shared.read_band(image_file, bands[b]['index'], sub=sub)
                data_shape = cdata_radiance.shape

                ## compute radiance
                cdata_radiance = cdata_radiance.astype(np.float32) * gain
                cdata_radiance += bias

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

            ## add current tile to outputs
            ofiles.append(ofile)

    return(ofiles, setu)
