## def l1_convert
## converts HICO L1 NC file to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-08-03
## modifications: 2021-12-31 (QV) new handling of settings
##                2022-01-04 (QV) added netcdf compression
##                2023-01-31 (QV) moved F0 import
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-16 (QV) use new gem NetCDF handling
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming

def l1_convert(inputfile, output = None, settings = None):
    import numpy as np
    import datetime, dateutil.parser, os
    import acolite as ac
    from netCDF4 import Dataset

    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    sensor = 'ISS_HICO'

    ## get sensor specific defaults
    setd = ac.acolite.settings.parse(sensor)
    ## set sensor default if user has not specified the setting
    for k in setd:
        if k not in ac.settings['user']: setu[k] = setd[k]
    ## end set sensor specific defaults

    verbosity = setu['verbosity']
    if output is None: output = setu['output']

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ofiles = []
    for file in inputfile:
        if output is None: output = os.path.dirname(file)
        hatts = ac.hico.attributes(file)

        ## get scene mid point time
        stime = dateutil.parser.parse('{}T{}'.format(hatts['Beginning_Date'], hatts['Beginning_Time']))
        etime = dateutil.parser.parse('{}T{}'.format(hatts['Ending_Date'], hatts['Ending_Time']))
        elapsed = (etime-stime).seconds
        time = stime + datetime.timedelta(seconds=elapsed/2)

        doy = int(time.strftime('%j'))
        d = ac.shared.distance_se(doy)

        gatts =  {}
        gatts['sensor'] = sensor
        gatts['isodate'] = time.isoformat()
        gatts['acolite_file_type'] = 'L1R'

        oname  = '{}_{}'.format(gatts['sensor'],  time.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## get F0 for radiance -> reflectance computation
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])

        gemo = ac.gem.gem(ofile, new = True)

        ## read geometry data
        ave = {}
        for ds in ['lat', 'lon', 'vza', 'vaa', 'sza', 'saa']:
            print('Reading HICO {}'.format(ds))
            data, att = ac.hico.read(file, ds)
            ave[ds] = np.nanmean(data)
            gemo.write(ds, data)

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

            if setu['output_lt']:
                ## write toa radiance
                gemo.write('Lt_{}'.format(bands[b]['wave_name']), cdata_radiance, ds_att = ds_att)
            ## write toa reflectance
            gemo.write('rhot_{}'.format(bands[b]['wave_name']), cdata, ds_att = ds_att)

        ## update gatts
        gemo.gatts = {k: gatts[k] for k in gatts}
        gemo.gatts_update()
        gemo.close()
        ofiles.append(ofile)
    return(ofiles, setu)
