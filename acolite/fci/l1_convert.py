## def l1_convert
## convert FCI nc bundle to ACOLITE L1R
## note that fcidecomp has to be installed, or that files have to be repacked using nccopy -F none
##
## written by Quinten Vanhellemont, RBINS
## 2025-05-13
## modifications: 2024-10-18 (QV) initial development
##                2025-05-12 (QV) added geolocation/geometry, radiance conversion and subsetting
##                2025-05-13 (QV) converted to l1_convert function

def l1_convert(inputfile, output = None, settings = None):
    import os, json
    import dateutil.parser, time
    import numpy as np
    import acolite as ac

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
            inputfile = inputfile.split(';') ## use semicolon as EUMETSAT uses commas in their file names
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ## VSWIR band naming from nc and rsr
    band_names = {'vis_04': 'VIS0.4', 'vis_05': 'VIS0.5',
                  'vis_06': 'VIS0.6_HR', 'vis_08': 'VIS0.8',
                  'vis_09': 'VIS0.9', 'nir_13': 'NIR1.3',
                  'nir_16': 'NIR1.6', 'nir_22': 'NIR2.2_HR'}

    ## run through inputfiles
    ofiles = []
    for fi, bundle in enumerate(inputfile):
        test = ac.fci.bundle_test(bundle)
        if test is None: continue

        ## retrieve fci files
        print('Loading FCI frames in repeat cycle for {}'.format(bundle))
        fci_files = ac.fci.find_files(bundle)

        ## list date_time_position files in the retrieved fci files
        fci_images = list(fci_files.keys())
        fci_images.sort()
        if len(fci_images) == 0:
            print('No FCI images found in {}'.format(bundle))
            continue

        ## run to image strips
        for dtp in fci_files.keys():
            ## determine subset once
            sub, lat, lon, vaa, vza = [None] * 5

            ## some hard coded values
            if 'HRFI' in dtp:
                datatype = 'HRFI'
                datasets = ['vis_06_hr', 'nir_22_hr']
                ## IR bands not used at the moment
                datasets_tir = ['ir_38_hr', 'ir_105_hr']
                full_dimensions = 22272, 22272
                ssd = 0.5
            elif 'FDHSI' in dtp:
                datatype = 'FDHSI'
                datasets = ['vis_04','vis_05','vis_06','vis_08','vis_09', 'nir_13','nir_16','nir_22']
                ## note that these may have different resolution, but not used at the moment
                datasets_tir = ['ir_38','wv_63','wv_73','ir_87','ir_97','ir_105','ir_123','ir_133']
                full_dimensions = 11136, 11136
                ssd = 1.0
            else:
                print('{} not recognised'.format(dtp))
                continue

            ## get attributes
            gatts = ac.shared.nc_gatts(fci_files[dtp][0][0])

            ## get date/time
            dt = dateutil.parser.parse(gatts['date_time_position'])

            ## get sensor type and position
            if 'nominalLongitude' in gatts:
                lon_0 = gatts['nominalLongitude']
            else:
                lon_0 = setu['fci_lon_0_default']
                print('Assuming default sub satellite longitude {} '.format(lon_0))

            platform, instrument = gatts['platform'], gatts['data_source']
            if platform in ['MTI1']:
                sensor = '{}_{}'.format('MTG-I{}'.format(platform[-1]), instrument)
            else:
                print('Platform {} not configured.'.format(platform))
                continue

            print(sensor, platform, instrument)

            ## get sensor specific defaults
            setd = ac.acolite.settings.parse(sensor)
            ## set sensor default if user has not specified the setting
            for k in setd:
                if k not in ac.settings['user']: setu[k] = setd[k]
            ## end set sensor specific defaults
            if output is None: output = setu['output']

            ## read rsr
            rsrd = ac.shared.rsr_dict(sensor)[sensor]

            #if not setu['use_provided_f0']:
            #    f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
            #    f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], rsrd['rsr'])

            ## global attributes
            gatts = {}
            gatts['sensor'] = sensor
            gatts['isodate'] = dt.isoformat()
            gatts['acolite_file_type'] = 'L1R'

            ## output name
            oname = '{}_{}'.format(gatts['sensor'], dt.strftime('%Y_%m_%d_%H_%M_%S'))
            if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
            ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
            #ofile_hrv = '{}/{}_{}_HRV.nc'.format(output, oname, gatts['acolite_file_type'])
            gatts['oname'] = oname
            gatts['ofile'] = ofile
            #gatts['ofile_hrv'] = ofile_hrv

            ## load geolocation and geometry
            if (lat is None) | (lon is None):
                lon, lat, vaa, vza = ac.seviri.geom(lon_0 = lon_0, ssd = ssd, instrument = instrument, geolocation = True, geometry = True)

            ## find image subset
            if setu['sub'] is not None:
                sub = setu['sub']
            elif setu['limit'] is not None:
                sub = ac.shared.geolocation_sub(lat, lon, setu['limit'])
                if sub is None:
                    print('Limit {} not in full disk.'.format(limit))
                    continue
            else:
                print('Warning: running without subsetting full disk image.')
                print('It is recommended using a subset sub=x, y, nx, ny or limit=S, W, N, E.')
                return

            ## subset geolocation
            column_range, row_range = None, None
            if sub is not None:
                column_range = sub[0], sub[0]+sub[2]
                row_range = sub[1], sub[1]+sub[3]

                ## subset geolocation and geometry
                lon = lon[row_range[0]:row_range[1], column_range[0]:column_range[1]]
                lat = lat[row_range[0]:row_range[1], column_range[0]:column_range[1]]
                vaa = vaa[row_range[0]:row_range[1], column_range[0]:column_range[1]]
                vza = vza[row_range[0]:row_range[1], column_range[0]:column_range[1]]

            ## get sun geometry
            print('Computing sun position for {}'.format(dt.isoformat()))
            spos = ac.shared.sun_position(dt, lon, lat)
            se_distance = spos['distance']
            sza, saa = spos['zenith'], spos['azimuth']
            del spos

            ## geom function already provides mask
            mask = np.isnan(lon)

            ## additional mask for sun below horizon
            ## can be made an option?
            mask = mask | (sza > 90)

            ## cosine sun zenith angle
            mus = np.cos(np.radians(sza))
            print('Cosine sun zenith angle shape: {}'.format(mus.shape))

            ## relative azimuth
            raa = np.abs(saa-vaa)
            raa[raa>180] = 360 - raa[raa>180]

            ## new gem handling
            gemo = ac.gem.gem(ofile, new = True)
            gemo.gatts = {k: gatts[k] for k in gatts}

            ## write position and angles
            gemo.write('lon', lon)
            del lon
            gemo.write('lat', lat)
            del lat
            gemo.write('vza', vza)
            del vza
            gemo.write('vaa', vaa)
            del vaa
            gemo.write('sza', sza)
            del sza
            gemo.write('saa', saa)
            del saa
            gemo.write('raa', raa)
            del raa

            ## run through bands
            for bi, ds in enumerate(band_names):
                band = band_names[ds]
                ds_att = {'wavelength':rsrd['wave_nm'][band]}

                ## read radiance
                print('Reading radiance for {} {}'.format(dtp, ds))
                data, irr = ac.fci.read_tiles(fci_files, dtp, ds, column_range = column_range, row_range = row_range)
                band_irr = irr * 1.0

                ### get F0
                #if setu['use_provided_f0']:
                #    band_irr = irr * 1.0
                #else:
                #    band_irr = f0d[band]
                ds_att['band_irr'] = band_irr
                print(band, band_irr)

                if setu['output_lt']:
                    ds = 'Lt_{}'.format(rsrd['wave_name'][band])
                    gemo.write(ds, data, ds_att = ds_att)

                ## convert to reflectance
                scale = (np.pi * se_distance ** 2) / (band_irr * mus)
                data *= scale

                ## mask negative TOA
                data[data <= 0] = np.nan
                data[mask] = np.nan

                ## write dataset
                ds = 'rhot_{}'.format(rsrd['wave_name'][band])
                gemo.write(ds, data, ds_att = ds_att)
                if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data.shape))
                data = None
            gemo.close()
            gemo = None
            if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
