## def l1_convert
## converts FORMOSAT5 bundle to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2022-04-12
## modifications: 2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-16 (QV) use new gem NetCDF handling
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming

def l1_convert(inputfile, output = None, settings = None):
    import numpy as np

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

    ofiles = []
    for bundle in inputfile:
        tiles, metafile = ac.formosat.bundle_test(bundle)
        meta = ac.formosat.metadata(metafile)
        if meta['sensor'] not in ['FORMOSAT5_RSI']: continue

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
        gatts['acolite_file_type'] = 'L1R'

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
        if setu['limit'] is None: sub = None

        ## run through tiles
        for ti, image_file in enumerate(tiles):
            oname  = '{}_{}'.format(gatts['sensor'],  dtime.strftime('%Y_%m_%d_%H_%M_%S'))
            if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
            ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
            gatts['oname'] = oname
            gatts['ofile'] = ofile

            gemo = ac.gem.gem(ofile, new = True)
            gemo.gatts = {k: gatts[k] for k in gatts}

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
                    gemo.write('lon', lon)
                    lon = None
                    if verbosity > 1: print('Wrote lon')
                    gemo.write('lat', lat)
                    lat = None
                    if verbosity > 1: print('Wrote lat')

            ## run through bands
            for bi, b in enumerate(band_names):
                if b == 'PAN': continue
                gain = meta['{}_PHYSICAL_GAIN'.format(b)]
                bias = meta['{}_PHYSICAL_BIAS'.format(b)]
                print('Computing rhot_{} for {}'.format(bands[b]['wave_name'], gatts['oname']))
                print(b, bi, gain, bias)

                ## read data
                cdata_radiance = ac.shared.read_band(image_file, bands[b]['index'], sub=sub)
                data_shape = cdata_radiance.shape

                ## compute radiance
                cdata_radiance = cdata_radiance.astype(np.float32) * gain
                cdata_radiance += bias

                if setu['output_lt']:
                    ## write toa radiance
                    gemo.write('Lt_{}'.format(bands[b]['wave_name']), cdata_radiance, ds_att = bands[b])

                ## compute reflectance
                cdata = cdata_radiance * (np.pi * gatts['se_distance'] * gatts['se_distance']) / (bands[b]['f0'] * mu0)
                cdata_radiance = None
                gemo.write('rhot_{}'.format(bands[b]['wave_name']), cdata, ds_att = bands[b])
                cdata = None
            gemo.close()
            ## add current tile to outputs
            ofiles.append(ofile)

    return(ofiles, setu)
