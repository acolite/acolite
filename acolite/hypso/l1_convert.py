## def l1_convert
## converts HYPSO file to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2023-03-09
## modifications: 2023-03-27 (QV) fixed issue with band naming

def l1_convert(inputfile, output = None, settings = {}, verbosity = 5):
    import numpy as np
    import datetime, dateutil.parser, os, copy
    import acolite as ac
    import h5py

    ## parse settings
    sensor = 'HYPSO1'
    setu = ac.acolite.settings.parse(sensor, settings=settings)
    verbosity = setu['verbosity']
    if output is None: output = setu['output']
    vname = setu['region_name']
    poly = setu['polygon']
    limit = setu['limit']

#     ## check if ROI polygon is given
#     clip, clip_mask = False, None
#     if poly is not None:
#         if os.path.exists(poly):
#             try:
#                 limit = ac.shared.polygon_limit(poly)
#                 print('Using limit from polygon envelope: {}'.format(limit))
#                 clip = True
#             except:
#                 print('Failed to import polygon {}'.format(poly))

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
        #gatts = ac.shared.nc_gatts(bundle)
        f = h5py.File(bundle, mode='r')

        ## get band information
        attributes = {k: f['/products']['Lt'].attrs[k] for k in f['/products']['Lt'].attrs.keys() if k != 'DIMENSION_LIST'}
        waves = attributes['wavelengths']
        fwhm = attributes['fwhm']

        rsr = ac.shared.rsr_hyper(waves, fwhm, step=0.1)
        rsrd = ac.shared.rsr_dict(rsrd={sensor:{'rsr':rsr}})

        ## create band dataset
        band_names = ['{}'.format(b) for b in range(0, len(waves))]
        #rsr = {'{}'.format(b): ac.shared.gauss_response(waves[bi], fwhm[bi], step=0.1) \
        #       for bi, b in enumerate(band_names)}
        #band_rsr= {b: {'wave': rsr[b][0]/1000, 'response': rsr[b][1]} for b in rsr}
        #f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], band_rsr)
        f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], rsr)
        bands = {}

        for bi, b in enumerate(band_names):
            #cwave = waves[bi]
            #swave = '{:.0f}'.format(cwave)
            cwave = rsrd[sensor]['wave_nm'][bi]
            swave = rsrd[sensor]['wave_name'][bi]
            bands[swave]= {'wave':cwave, 'wavelength':cwave, 'wave_mu':cwave/1000.,
                           'wave_name':swave, 'width': fwhm[bi],
                           'rsr': rsr[bi],'f0': f0d[bi]}

        ## time
        utime = f['/navigation/']['unixtime'][:]

        ## geometry
        vza = f['/navigation/']['sensor_zenith'][:]
        vaa = f['/navigation/']['sensor_azimuth'][:]

        sza = f['/navigation/']['solar_zenith'][:]
        saa = f['/navigation/']['solar_azimuth'][:]

        ## compute relative azimuth angle
        raa = np.abs(saa-vaa)
        ## raa along 180 degree symmetry
        tmp = np.where(raa>180)
        raa[tmp]=np.abs(raa[tmp] - 360)

        ## read lat and lon
        lat = f['/navigation/']['latitude'][:]
        lon = f['/navigation/']['longitude'][:]

        ## metadata
        start_time = datetime.datetime.utcfromtimestamp(utime[0])
        stop_time = datetime.datetime.utcfromtimestamp(utime[-1])
        #start_time = dateutil.parser.parse(t0)
        #stop_time = dateutil.parser.parse(t1)

        ## date time
        dt = start_time + (stop_time-start_time)
        isodate = dt.isoformat()

        ## compute scene center sun position
        clon = np.nanmedian(lon)
        clat = np.nanmedian(lat)
        cspos = ac.shared.sun_position(dt, clon, clat)
        #cossza = np.cos(np.radians(spos['zenith']))

        spos = ac.shared.sun_position(dt, lon, lat)
        d = spos['distance']
        #sza = spos['zenith']
        #saa = spos['azimuth']
        spos = None

        ## output attributes
        gatts = {}
        gatts['sensor'] = sensor
        gatts['acolite_file_type'] = 'L1R'
        gatts['isodate'] = dt.isoformat()

        ## centre sun position
        gatts['doy'] = dt.strftime('%j')
        gatts['se_distance'] = cspos['distance']

        gatts['sza'] = np.nanmean(sza)
        gatts['saa'] = np.nanmean(saa)
        gatts['mus'] = np.cos(np.radians(gatts['sza']))
        gatts['vza'] = np.nanmean(vza)
        gatts['vaa'] = np.nanmean(vaa)

        gatts['raa'] = np.abs(gatts['saa']-gatts['vaa'])
        while gatts['raa'] > 180: gatts['raa'] = np.abs(gatts['raa']-360)

        ## store band info in gatts
        gatts['band_waves'] = waves # [bands[w]['wave'] for w in bands]
        gatts['band_widths'] = fwhm # [bands[w]['width'] for w in bands]

        ## output file
        obase  = '{}_{}_{}'.format(gatts['sensor'],  dt.strftime('%Y_%m_%d_%H_%M_%S'), gatts['acolite_file_type'])
        if not os.path.exists(output): os.makedirs(output)
        ofile = '{}/{}.nc'.format(output, obase)
        gatts['obase'] = obase

        new = True
        nc_projection = None
        ## write lat/lon
        if (setu['output_geolocation']) & (new):
            if verbosity > 1: print('Writing geolocation lon/lat')
            ac.output.nc_write(ofile, 'lon', lon, attributes=gatts, new=new,
                               double=True, nc_projection=nc_projection,
                               netcdf_compression=setu['netcdf_compression'],
                               netcdf_compression_level=setu['netcdf_compression_level'])
            if verbosity > 1: print('Wrote lon ({})'.format(lon.shape))

            ac.output.nc_write(ofile, 'lat', lat, double=True,
                                netcdf_compression=setu['netcdf_compression'],
                                netcdf_compression_level=setu['netcdf_compression_level'])
            if verbosity > 1: print('Wrote lat ({})'.format(lat.shape))
            new=False

        ## write geometry
        if (setu['output_geometry']):
            if verbosity > 1: print('Writing geometry')
            ac.output.nc_write(ofile, 'vza', vza, attributes=gatts, new=new,
                               double=True, nc_projection=nc_projection,
                               netcdf_compression=setu['netcdf_compression'],
                               netcdf_compression_level=setu['netcdf_compression_level'])
            if verbosity > 1: print('Wrote vza ({})'.format(vza.shape))
            vza = None
            new=False

            ac.output.nc_write(ofile, 'vaa', vaa, attributes=gatts, new=new,
                               double=True, nc_projection=nc_projection,
                               netcdf_compression=setu['netcdf_compression'],
                               netcdf_compression_level=setu['netcdf_compression_level'])
            if verbosity > 1: print('Wrote vaa ({})'.format(vaa.shape))
            vaa = None

            ac.output.nc_write(ofile, 'sza', sza, attributes=gatts, new=new,
                               double=True, nc_projection=nc_projection,
                               netcdf_compression=setu['netcdf_compression'],
                               netcdf_compression_level=setu['netcdf_compression_level'])
            if verbosity > 1: print('Wrote sza ({})'.format(sza.shape))

            ac.output.nc_write(ofile, 'saa', saa, attributes=gatts, new=new,
                               double=True, nc_projection=nc_projection,
                               netcdf_compression=setu['netcdf_compression'],
                               netcdf_compression_level=setu['netcdf_compression_level'])
            if verbosity > 1: print('Wrote saa ({})'.format(saa.shape))
            saa = None

            ac.output.nc_write(ofile, 'raa', raa, attributes=gatts, new=new,
                               double=True, nc_projection=nc_projection,
                               netcdf_compression=setu['netcdf_compression'],
                               netcdf_compression_level=setu['netcdf_compression_level'])
            if verbosity > 1: print('Wrote raa ({})'.format(raa.shape))
            raa = None

        ## cosine of sun angle
        cossza = np.cos(np.radians(sza))
        sza = None

        ## read data cube
        data = f['/products']['Lt'][:]

        ## close file
        f.close()
        f = None

        nbands = data.shape[2]
        data_dimensions = data.shape[0], data.shape[1]

        ## run through bands and store rhot
        for bi, b in enumerate(bands):
            ## get dataset attributes
            ds_att = {k:bands[b][k] for k in bands[b] if k not in ['rsr']}

            ## copy radiance
            cdata_radiance = data[:,:,bi]

            if setu['output_lt']:
                ## write toa radiance
                ac.output.nc_write(ofile, 'Lt_{}'.format(bands[b]['wave_name']), cdata_radiance,
                                      dataset_attributes = ds_att, new = new, attributes = gatts,
                                      netcdf_compression=setu['netcdf_compression'],
                                      netcdf_compression_level=setu['netcdf_compression_level'],
                                      netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
                new = False
                print('Wrote Lt_{}'.format(bands[b]['wave_name']))

            ## compute reflectance
            cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * cossza)
            cdata_radiance = None

            ## write toa reflectance
            ac.output.nc_write(ofile, 'rhot_{}'.format(bands[b]['wave_name']), cdata,
                                dataset_attributes = ds_att, new = new, attributes = gatts,
                                netcdf_compression=setu['netcdf_compression'],
                                netcdf_compression_level=setu['netcdf_compression_level'],
                                netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
            cdata = None
            new = False
            print('Wrote rhot_{}'.format(bands[b]['wave_name']))
        data = None
        ofiles.append(ofile)
    return(ofiles, setu)
