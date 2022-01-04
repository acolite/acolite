## def l1_convert
## converts HYPERION L1T bundle to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-08-04
## modifications:  2021-12-31 (QV) new handling of settings
##                2022-01-04 (QV) added netcdf compression

def l1_convert(inputfile, output = None, settings = {}, verbosity=5):
    import numpy as np
    import datetime, dateutil.parser, os, copy
    import acolite as ac
    from netCDF4 import Dataset

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

    ## get HYPERION wavelengths
    hypf = ac.config['data_dir'] + '/HYPERION/Hyperion_cen_fwhm.dat'
    hypd = np.loadtxt(hypf, encoding='utf-8')
    nbands = hypd.shape[0]

    ofiles = []
    for file in inputfile:
        ## find metadata
        metadata = ac.hyperion.metadata(file)
        if metadata['PRODUCT_METADATA']['SENSOR_ID'] != 'HYPERION': continue

        ## output attributes
        gatts = {}

        ## get dimensions
        nrow=int(metadata['PRODUCT_METADATA']['PRODUCT_LINES'])
        ncol=int(metadata['PRODUCT_METADATA']['PRODUCT_SAMPLES'])
        dim = (nrow,ncol)

        ## scene start and end time, time per row
        stime = dateutil.parser.parse('{}T{}'.format(metadata['PRODUCT_METADATA']['ACQUISITION_DATE'],
                                                     metadata['PRODUCT_METADATA']['START_TIME'].split()[-1]))
        etime = dateutil.parser.parse('{}T{}'.format(metadata['PRODUCT_METADATA']['ACQUISITION_DATE'],
                                                     metadata['PRODUCT_METADATA']['END_TIME'].split()[-1]))
        tdiff = (etime-stime).seconds
        trow = tdiff/nrow
        time = stime + datetime.timedelta(seconds = tdiff/2)

        doy = int(time.strftime('%j'))
        d = ac.shared.distance_se(doy)
        gatts['isodate'] = time.isoformat()

        ## more metadata
        sensor = metadata['PRODUCT_METADATA']['SENSOR_ID']
        satellite = metadata['PRODUCT_METADATA']['SPACECRAFT_ID']

        gatts['sensor'] = '{}_{}'.format(satellite, sensor)

        ## get sensor specific settings
        setu = ac.acolite.settings.parse(gatts['sensor'], settings=settings)
        verbosity = setu['verbosity']
        limit = setu['limit']
        vname = setu['region_name']
        output_lt=setu['output_lt']
        if output is None: output = setu['output']

        gatts['sza'] = 90-float(metadata["PRODUCT_PARAMETERS"]['SUN_ELEVATION'])
        gatts['vza'] = np.abs(float(metadata["PRODUCT_PARAMETERS"]['SENSOR_LOOK_ANGLE']))
        if 'raa' not in gatts: gatts['raa'] = 0
        mu0 = np.cos(gatts['sza']*(np.pi/180))
        muv = np.cos(gatts['vza']*(np.pi/180))

        cell_size = float(metadata['PROJECTION_PARAMETERS']['GRID_CELL_SIZE'])
        projection = metadata['PROJECTION_PARAMETERS']['MAP_PROJECTION']
        datum = metadata['PROJECTION_PARAMETERS']['REFERENCE_DATUM']
        ellipsoid = metadata['PROJECTION_PARAMETERS']['REFERENCE_ELLIPSOID']
        zone = metadata['UTM_PARAMETERS']['ZONE_NUMBER']

        ## get band indices and scaling
        bands_swir = [i for i in range(71,243)]
        bands_vnir = [i for i in range(1,71)]
        scaling_swir = float(metadata['RADIANCE_SCALING']['SCALING_FACTOR_SWIR']) ## bands 71-242 (W/(m^2 sr um)
        scaling_vnir = float(metadata['RADIANCE_SCALING']['SCALING_FACTOR_VNIR']) ## bands 1-70 (W/(m^2 sr um)

        if output is None:
            odir = os.path.dirname(file)
        else:
            odir = output

        obase  = '{}_{}_L1R'.format(gatts['sensor'],  time.strftime('%Y_%m_%d_%H_%M_%S'))
        if not os.path.exists(odir): os.makedirs(odir)
        ofile = '{}/{}.nc'.format(odir, obase)
        gatts['obase'] = obase

        ## add band info
        gatts['band_waves'] = hypd[:, 1]
        gatts['band_widths'] = hypd[:, 2]

        ## make rsr and bands dataset
        rsr = {'{}'.format(b): ac.shared.gauss_response(gatts['band_waves'][b], gatts['band_widths'][b], step=0.1)
               for b in range(len(gatts['band_waves']))}
        band_rsr = {b: {'wave': rsr[b][0]/1000, 'response': rsr[b][1]}  for b in rsr}
        f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], band_rsr)

        ## make bands dataset
        bands = {}
        for bi, b in enumerate(band_rsr):
            cwave = gatts['band_waves'][bi]
            swave = '{:.0f}'.format(cwave)
            bands[b]= {'wave':cwave, 'wavelength':cwave, 'wave_mu':cwave/1000.,
                       'wave_name':'{:.0f}'.format(cwave),
                       'width': gatts['band_widths'][bi],
                       'rsr': band_rsr[b],'f0': f0d[b]}

        ## generate new file
        new = True

        ## get scene projection info
        dct = ac.hyperion.projection(metadata)

        ## find region of interest
        sub = None
        if (sub is None) & (limit is not None):
            dct_sub = ac.shared.projection_sub(dct, limit, four_corners=True)
            if dct_sub['out_lon']:
                if verbosity > 1: print('Longitude limits outside {}'.format(file))
                continue
            if dct_sub['out_lat']:
                if verbosity > 1: print('Latitude limits outside {}'.format(file))
                continue
            sub = dct_sub['sub']
            dct_prj = {k:dct_sub[k] for k in dct_sub}
        else:
            dct_prj = {k:dct[k] for k in dct}

        ## compute longitude and latitude
        lon, lat = ac.shared.projection_geo(dct_prj)
        ac.output.nc_write(ofile, 'lon', lon, attributes=gatts, new=new, double=True,
                            netcdf_compression=setu['netcdf_compression'],
                            netcdf_compression_level=setu['netcdf_compression_level'])
        if verbosity > 1: print('Wrote lon')
        lon = None

        ac.output.nc_write(ofile, 'lat', lat, double=True,
                            netcdf_compression=setu['netcdf_compression'],
                            netcdf_compression_level=setu['netcdf_compression_level'])
        if verbosity > 1: print('Wrote lat')
        lat = None
        new = False

        for b,band in enumerate(bands_vnir+bands_swir):
            if band in bands_vnir: det = 'VNIR'
            if band in bands_swir: det = 'SWIR'

            btag = '{}'.format(band-1)

            ## output datasets
            ds_att = {k:bands[btag][k] for k in bands[btag] if k not in ['rsr']}

            print('Processing HYPERION band {} on {} detector rhot_{}'.format(b+1, det, ds_att['wave_name']))

            ## read DN and convert to toa radiance
            band_file = metadata['PRODUCT_METADATA']['BAND{}_FILE_NAME'.format(band)]
            cdata_radiance = ac.shared.read_band('{}/{}'.format(file, band_file), sub = sub).astype(np.float32)
            if det == 'VNIR': cdata_radiance /= scaling_vnir
            if det == 'SWIR': cdata_radiance /= scaling_swir
            cdata_radiance[cdata_radiance == 0] = np.nan

            ## write toa radiance
            if output_lt:
                ac.output.nc_write(ofile, 'Lt_{}'.format(ds_att['wave_name']), cdata_radiance,
                                   attributes = gatts, dataset_attributes = ds_att, new=new,
                                   netcdf_compression=setu['netcdf_compression'],
                                   netcdf_compression_level=setu['netcdf_compression_level'],
                                   netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
                new = False

            ## compute toa reflectance
            cdata = (np.pi*cdata_radiance*d*d) / (ds_att['f0'] * mu0)

            ## write toa reflectance
            ac.output.nc_write(ofile, 'rhot_{}'.format(ds_att['wave_name']), cdata,
                               attributes = gatts, dataset_attributes = ds_att, new=new,
                               netcdf_compression=setu['netcdf_compression'],
                               netcdf_compression_level=setu['netcdf_compression_level'],
                               netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
            new = False



        ofiles.append(ofile)
    return(ofiles, setu)
