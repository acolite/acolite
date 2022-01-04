## def l1_convert
## converts HICO L1 NC file to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-08-03
## modifications: 2021-12-31 (QV) new handling of settings
##                2022-01-04 (QV) added netcdf compression

def l1_convert(inputfile, output = None, settings = {}, verbosity=5):
    import numpy as np
    import datetime, dateutil.parser, os
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

    ofiles = []
    for file in inputfile:
        hatts = ac.hico.attributes(file)

        ## get scene mid point time
        stime = dateutil.parser.parse('{}T{}'.format(hatts['Beginning_Date'], hatts['Beginning_Time']))
        etime = dateutil.parser.parse('{}T{}'.format(hatts['Ending_Date'], hatts['Ending_Time']))
        elapsed = (etime-stime).seconds
        time = stime + datetime.timedelta(seconds=elapsed/2)

        doy = int(time.strftime('%j'))
        d = ac.shared.distance_se(doy)

        if output is None:
            odir = os.path.dirname(file)
        else:
            odir = output

        gatts =  {}
        gatts['sensor'] = 'ISS_HICO'
        gatts['isodate'] = time.isoformat()

        ## get sensor specific settings
        setu = ac.acolite.settings.parse(gatts['sensor'], settings=settings)
        verbosity = setu['verbosity']
        vname = setu['region_name']
        output_lt=setu['output_lt']
        if output is None: output = setu['output']

        obase  = '{}_{}_L1R'.format(gatts['sensor'],  time.strftime('%Y_%m_%d_%H_%M_%S'))
        if not os.path.exists(odir): os.makedirs(odir)
        ofile = '{}/{}.nc'.format(odir, obase)
        gatts['obase'] = obase

        ## read geometry data
        new = True
        ave = {}
        for ds in ['lat', 'lon', 'vza', 'vaa', 'sza', 'saa']:
            print('Reading HICO {}'.format(ds))
            data, att = ac.hico.read(file, ds)
            ave[ds] = np.nanmean(data)
            ac.output.nc_write(ofile, ds, data, new=new, attributes=gatts,
                                netcdf_compression=setu['netcdf_compression'],
                                netcdf_compression_level=setu['netcdf_compression_level'])
            new = False

        if 'sza' not in gatts: gatts['sza'] = ave['sza']
        if 'vza' not in gatts: gatts['vza'] = ave['vza']
        if 'saa' not in gatts: gatts['saa'] = ave['saa']
        if 'vaa' not in gatts: gatts['vaa'] = ave['vaa']

        if 'raa' not in gatts:
            raa_ave = abs(gatts['saa'] - gatts['vaa'])
            while raa_ave >= 180: raa_ave = abs(raa_ave-360)
            gatts['raa'] = raa_ave

        mu0 = np.cos(gatts['sza']*(np.pi/180))
        muv = np.cos(gatts['vza']*(np.pi/180))

        ## read hico Lt
        data, att = ac.hico.read(file, 'lt')

        gatts['band_waves'] = att['wavelengths']
        gatts['band_widths'] = att['fwhm']

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
                           'rsr': band_rsr[b],
                           'f0': f0d[b]}

        for bi, b in enumerate(band_rsr):
            print('Reading HICO rhot_{}'.format(bands[b]['wave_name']))
            cdata_radiance = data[:,:,bi]
            cdata = cdata_radiance * (np.pi * d * d) / (bands[b]['f0'] * mu0)

            ## output datasets
            ds_att = {k:bands[b][k] for k in bands[b] if k not in ['rsr']}

            if output_lt:
                ## write toa radiance
                ac.output.nc_write(ofile, 'Lt_{}'.format(bands[b]['wave_name']), cdata_radiance, dataset_attributes = ds_att,
                            netcdf_compression=setu['netcdf_compression'],
                            netcdf_compression_level=setu['netcdf_compression_level'],
                            netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])

            ## write toa reflectance
            ac.output.nc_write(ofile, 'rhot_{}'.format(bands[b]['wave_name']), cdata, dataset_attributes = ds_att,
                            netcdf_compression=setu['netcdf_compression'],
                            netcdf_compression_level=setu['netcdf_compression_level'],
                            netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])

        ## update gatts
        with Dataset(ofile, 'a') as nc:
            for key in gatts.keys():
                if gatts[key] is not None:
                    try:
                        nc.setncattr(key, gatts[key])
                    except:
                        print('Failed to write attribute: {}'.format(key))

        ofiles.append(ofile)
    return(ofiles, setu)
