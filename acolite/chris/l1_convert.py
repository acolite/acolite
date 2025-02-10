## def l1_convert
## converts CHRIS HDF file to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-06-08
## modifications: 2021-12-31 (QV) new handling of settings
##                2022-01-04 (QV) added netcdf compression
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-16 (QV) use new gem NetCDF handling
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming

def l1_convert(inputfile, output = None, settings = None):
    from pyhdf.SD import SD,SDC
    import datetime
    import numpy as np
    import acolite as ac

    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    sensor = 'CHRIS'

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

    ## get CHRIS interband calibration
    if setu['chris_interband_calibration']: ibcal = ac.chris.interband_calibration()

    ofiles = []
    for bundle in inputfile:
        if output is None: output = os.path.dirname(bundle)

        ## read gains and get bands
        gains, mode_info = ac.chris.vdata(bundle)
        nbands = len(mode_info)
        waves = [mode_info[b]['WlMid'] for b in range(0,nbands)]
        widths = [mode_info[b]['BWidth'] for b in range(0,nbands)]

        ## open HDF
        hdf = SD(bundle, SDC.READ)
        datasets = hdf.datasets()
        attributes = hdf.attributes()

        ## extract attributes
        sensor = attributes['Sensor Type']
        year, month, day = attributes['Image Date'].split('-')
        hour, minutes, seconds = attributes['Calculated Image Centre Time'].split(':')

        ## date/time
        tc = datetime.datetime(int(year), int(month), int(day), int(hour), int(minutes), int(seconds))
        isodate = tc.isoformat()
        doy = int(tc.strftime('%j'))

        ## target lat/lon
        tlat = float(attributes['Target Latitude'])
        tlon = float(attributes['Target Longitude'])

        ## compute sun position
        sd = ac.shared.sun_position(tc, tlon, tlat)
        sza_ = sd['zenith']
        saa_ = sd['azimuth']
        se_distance = sd['distance']

        #if 'Solar Azimuth Angle' in attributes:
        #    saa = float(attributes['Solar Azimuth Angle'])
        #else:
        #    saa_ = saa_[0]

        ## do angle calculations
        mza = float(attributes['Minimum Zenith Angle'])
        flyby_za = float(attributes['Fly-by Zenith Angle'])
        Palt = float(attributes['Platform Altitude'])
        Talt = float(attributes['Target Altitude'])
        Tlat = float(attributes['Target Latitude'])
        vza_, vaa_ = ac.chris.view_geometry(mza, flyby_za, Palt, Talt, Tlat)

        ## get view geometry
        vza = abs(vza_) if attributes['Observation Zenith Angle'] == 'unknown' \
                else float(attributes['Observation Zenith Angle'])
        if vza == 0.0: vza = 0.001
        vaa = abs(vaa_) if attributes['Observation Azimuth Angle'] == 'unknown' \
                else float(attributes['Observation Azimuth Angle'])
        muv = np.cos(vza*(np.pi/180))

        ## get sun geometry
        saa = abs(saa_[0]) if 'Solar Azimuth Angle' not in attributes \
                else float(attributes['Solar Azimuth Angle'])
        sza = float(attributes['Solar Zenith Angle'])
        mu0 = np.cos(sza*(np.pi/180))

        ## calculate relative azimuth angle
        raa = np.abs(saa-vaa)
        if raa > 180: raa = np.abs(360-raa)

        ## get bands and mode
        nbands = int(attributes['Number of Bands'])
        ncol = int(attributes['Number of Samples'])
        nrow = int(attributes['Number of Ground Lines'])
        mode = int(attributes['CHRIS Mode'])
        temp = float(attributes['CHRIS Temperature'])
        satellite = 'PROBA1'

        gatts = {'sensor':sensor, 'satellite':satellite, 'inputfile': bundle, \
                 'isodate':isodate, 'acolite_file_type': 'L1R',
                 'band_waves': waves, 'band_widths': widths,
                 'sza': sza, 'vza': vza, 'raa': raa, 'vaa': vaa, 'saa': saa,
                 'lat': tlat, 'lon': tlon, 'FlyByZenith':flyby_za}

        oname  = '{}_{}_{}_FlyByZenith_{}{}'.format(satellite, sensor,
                                                 tc.strftime('%Y_%m_%d_%H_%M_%S'),
                                                 '-' if gatts['FlyByZenith'] < 0 else '+',
                                                 '{:.0f}'.format(abs(gatts['FlyByZenith'])).zfill(2))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## make RSR
        rsr = ac.shared.rsr_hyper(gatts['band_waves'], gatts['band_widths'])

        #rsr = {}
        #for b in range(nbands):
        #    wave, resp = ac.shared.gauss_response(waves[b], widths[b], step=0.25)
        #    rsr[b] = {'wave':wave/1000, 'response': resp}

        ## get F0 per band
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        f0_bands = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], rsr)

        ## determine mode
        if mode==1:
            mode_name = 'Full swath width, 62 spectral bands, 773nm / 1036nm, nadir ground sampling distance 34m @ 556km'
        if mode==2:
            mode_name = 'WATER BANDS: Full swath width, 18 spectral bands, nadir ground sampling distance 17m @ 556km'
        if mode==3:
            mode_name = 'LAND CHANNELS: Full swath width, 18 spectral bands, nadir ground sampling distance 17m @ 556km'
        if mode==4:
            mode_name = 'CHLOROPHYL BAND SET: Full swath width, 18 spectral bands, nadir ground sampling distance 17m @ 556km'
        if mode==5:
            mode_name = 'LAND CHANNELS: Half swath width, 37 spectral bands, nadir ground sampling distance 17m @ 556km'
        if mode == 1:
            overscan_s = [0,2]
            dark_s = [2,2]
            blank_s = [4,2]
            image = [6,374]
            dark_e = [380,2]
            overscan_e = [382,2]
            padding = [384,382]
        elif mode in [2,3,4]:
            overscan_s = [0,4]
            dark_s = [4,4]
            blank_s = [8,4]
            image = [12,744]
            dark_e = [756,4]
            overscan_e = [760,6]
            padding = [0,0]
        elif mode == 5:
            overscan_s = [0,4]
            dark_s = [4,4]
            blank_s = [8,4]
            image = [12,370]
            overscan_e = [382,384]
            padding = [0,0]
        masks = {0:'Useful', 1:'Reset', 2:'Saturated'}

        ## select interband calibration
        cal = None
        if setu['chris_interband_calibration']:
            for c in ibcal:
                if (tc >= ibcal[c]['range'][0]) & (tc <= ibcal[c]['range'][1]):
                    cal = ibcal[c]['data']

        ## write TOA reflectance
        if verbosity > 1: print('Converting bands')

        ## set up output file
        gemo = ac.gem.gem(ofile, new=True)
        gemo.gatts = {k: gatts[k] for k in gatts}

        for b in range(nbands):
            wave_f = mode_info[b]['WlMid']
            width = mode_info[b]['BWidth']
            wave = '{:.0f}'.format(wave_f)

            if mode_info[b]['Gain'] in gains:
                gain_f = gains[mode_info[b]['Gain']]
            else:
                gain_f = 1

            for dataset in datasets:
                data = hdf.select(dataset)
                if dataset == 'RCI Image':
                    cur_data = data[b,0:nrow,image[0]:image[0]+image[1]] * 1.0
                if dataset == 'Saturation/Reset Mask':
                    cur_mask = data[b,0:nrow,image[0]:image[0]+image[1]]
            cur_data[cur_mask != 0] = np.nan


            ds_att = {'chris_gain': gain_f,
                      'wavelength': wave_f,
                      'width': width,
                      #'rsr_wave': rsr[b]['wave'],
                      #'rsr_response': rsr[b]['response']
                      'band_index': b, 'f0': f0_bands[b]}

            ## do interband calibration
            if (setu['chris_interband_calibration']) & (cal is not None):
                btag = 'band_{}'.format(b+1)
                if btag in cal:
                    bcal = cal[btag]['cal']
                    cur_data *= bcal
                    ds_att['interband_calibration'] = bcal
                    print(wave_f, bcal)

            ## output TOA radiance
            if setu['output_lt']:
                ds = 'L_{}'.format(wave)
                gemo.write(ds, cur_data, ds_att = ds_att)
                if verbosity > 2: print('Wrote {} to {}'.format(ds, ofile))

            ## output TOA reflectance
            cur_data *= (np.pi*se_distance*se_distance) / (f0_bands[b]*1000 * mu0)
            ds = 'rhot_{}'.format(wave)
            gemo.write(ds, cur_data, ds_att = ds_att)
            if verbosity > 2: print('Wrote {} to {}'.format(ds, ofile))
        hdf = None
        gemo.close()

        ## apply noise reduction
        if setu['chris_noise_reduction']:
            print('Applying CHRIS Noise Reduction')
            ofile = ac.chris.noise_reduction(ofile)

        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
