## def l1_convert
## convert Himawari HSD bundle to ACOLITE L1R
##
## written by Quinten Vanhellemont, RBINS
## 2025-05-19
## modifications: 2025-05-20 (QV) mask and resample at segment level, added crop at segment level
##                2025-05-22 (QV) only delete lon, lat, vaa, vza if last scene for the sensor
##                                get observation time from (sub)scene centre

def l1_convert(inputfile, output = None, settings = None):
    import os, json
    import datetime, dateutil.parser, time
    import scipy.ndimage
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

    ## run through inputfiles
    ofiles = []
    for fi, bundle in enumerate(inputfile):
        ## identify Himawari
        fd = ac.himawari.bundle_test(bundle)
        if len(fd) == 0: continue

        ## run through the files we found
        for platform in fd:
            ## import rsr
            sensor = '{}_AHI'.format(platform)
            rsrd = ac.shared.rsr_dict(sensor)[sensor]

            ## get sensor specific defaults
            setd = ac.acolite.settings.parse(sensor)
            ## set sensor default if user has not specified the setting
            for k in setd:
                if k not in ac.settings['user']: setu[k] = setd[k]
            ## end set sensor specific defaults
            if output is None: output = setu['output']

            ## target resolution
            ssd = setu['ahi_target_res']
            lon_0 = setu['ahi_lon_0_default']

            if ssd == 0.5:
                shape = 22000
            elif ssd == 1:
                shape = 11000
            elif ssd == 2:
                shape = 5500
            else:
                print('ahi_target_res={} not supported'.format(setu['ahi_target_res']))
                continue

            ## load geolocation and geometry
            #lon, lat, vaa, vza = ac.seviri.geom(lon_0 = lon_0, ssd = ssd, instrument = 'AHI', geolocation = True, geometry = True)
            lon, lat = ac.seviri.geom(lon_0 = lon_0, ssd = ssd, instrument = 'AHI', geolocation = True, geometry = False)
            geom_shape = lon.shape

            ## find image subset
            sub = None
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
                #return

            ## get F0 for radiance -> reflectance computation
            f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
            f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], rsrd['rsr'])

            ## subset geolocation
            column_range, row_range = None, None
            if sub is not None:
                column_range = sub[0], sub[0]+sub[2]
                row_range = sub[1], sub[1]+sub[3]
                ## subset geolocation
                lon = lon[row_range[0]:row_range[1], column_range[0]:column_range[1]]
                lat = lat[row_range[0]:row_range[1], column_range[0]:column_range[1]]

            ## load angles
            vaa, vza = ac.seviri.geom(lon_0 = lon_0, ssd = ssd, instrument = 'AHI', geolocation = False, geometry = True, sub = sub)

            ## run through dates for current sensor
            for date_idx, date in enumerate(fd[platform]):

                gatts = {}

                ## find location
                for bi, band in enumerate(rsrd['rsr_bands']):
                    band_name = 'B{}'.format(band.zfill(2))
                    if band_name not in fd[platform][date]:
                        print('Band {} missing'.format(band_name))
                        continue

                    ## output dataset attributes
                    ds_att = {'wavelength':rsrd['wave_nm'][band], 'f0': f0d[band]}

                    ## find segments
                    seg_i = [int(di[1:3]) for fi, di in enumerate(fd[platform][date][band_name])]
                    seg_n = [int(di[3:]) for fi, di in enumerate(fd[platform][date][band_name])]

                    if len(seg_i) != seg_n[0]:
                        print('Not enough segments for band {}: {}'.format(band_name, ','.join([str(v) for v in seg_i])))
                        continue

                    print('Segments for band {}: {}'.format(band_name, ','.join([str(v) for v in seg_i])))
                    segments = list(fd[platform][date][band_name].keys())
                    segments.sort()

                    print(band_name)

                    ## load band data
                    band_data = []
                    mask_data = []
                    ## track line times
                    observation_lines = []
                    observation_times = []
                    for si, segment in enumerate(segments):
                        segment_file = fd[platform][date][band_name][segment]['path']

                        ## read header
                        header = ac.himawari.segment_parse(segment_file, parse_data = False)

                        ## get observation times
                        for i in range(header[9]['number_of_observation_times']):
                            oi = i + 1
                            observation_lines.append(header[9]['line_number_{}'.format(oi)])
                            observation_times.append(header[9]['observation_time_{}'.format(oi)])

                        ## do subsetting
                        if sub is not None:
                            #print(band_name, header[2]['number_of_columns'], shape)
                            ## add one offset to start at 0
                            if header[2]['number_of_columns'] != shape:
                                factor = header[2]['number_of_columns']/shape
                                l0 = int((header[7]['first_line_number_of_image_sequence'] - 1) / factor)
                                l1 = int((l0 + (header[2]['number_of_lines'])/ factor))
                            else:
                                factor = 1
                                l0 = header[7]['first_line_number_of_image_sequence']-1
                                l1 = l0 + header[2]['number_of_lines']-1

                            row_range_segment = [row_range[0]-l0, row_range[1]-l0]
                            if row_range_segment[0] < 0: row_range_segment[0] = 0
                            if row_range_segment[1] > l1: row_range_segment[1] = l1

                            if (row_range_segment[0] > l1) | (row_range_segment[1] < 0):
                                print('Skipping segment {} for subset {}:{} {}:{}.'.format(segment, row_range[0], row_range[1], column_range[0], column_range[1]))
                                continue

                            ## read data and compute mask
                            data = ac.himawari.segment_parse(segment_file, header=header)
                            mask = data == header[5]['count_value_of_pixels_outside_scan_area']
                            mask = mask | data == header[5]['count_value_of_error_pixels']

                            ## resample using segment shape
                            if header[2]['number_of_columns'] != shape:
                                data = scipy.ndimage.zoom(data, zoom=1/factor, order=1)
                                mask = scipy.ndimage.zoom(mask, zoom=1/factor, order=1)

                            ## crop to subset in current segment
                            data = data[row_range_segment[0]:row_range_segment[1], column_range[0]:column_range[1]]
                            mask = mask[row_range_segment[0]:row_range_segment[1], column_range[0]:column_range[1]]

                        else:
                            ## read data and compute mask
                            data = ac.himawari.segment_parse(segment_file, header=header)
                            mask = band_data == header[5]['count_value_of_pixels_outside_scan_area']
                            mask = mask | band_data == header[5]['count_value_of_error_pixels']

                        ## append data
                        band_data.append(data)
                        mask_data.append(mask)
                        del data, mask

                    ## set up gatts and create new file
                    if len(gatts) == 0:
                        ## get time
                        ## MJD starts Nov. 17 1858
                        dt0 = datetime.datetime(1858, 11, 17, 0, 0, 0)

                        ## get observation time
                        observation_lines = np.asarray(observation_lines)
                        observation_times = np.asarray(observation_times)
                        factor = header[2]['number_of_columns']/shape
                        if sub is not None: centre_pixel = (row_range[1]-row_range[0])/2
                        else: centre_pixel = shape/2
                        centre_time = np.interp(centre_pixel, observation_lines, observation_times)
                        dt = dt0 + datetime.timedelta(days=centre_time)

                        ## get up gatts
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

                        ## get sun geometry
                        print('Computing sun position for {}'.format(dt.isoformat()))
                        spos = ac.shared.sun_position(dt, lon, lat)
                        se_distance = spos['distance']
                        sza, saa = spos['zenith'], spos['azimuth']
                        del spos

                        sza_mask = sza > 90

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
                        if date_idx == len(fd[platform])-1: del lon
                        gemo.write('lat', lat)
                        if date_idx == len(fd[platform])-1: del lat
                        gemo.write('vza', vza)
                        if date_idx == len(fd[platform])-1: del vza
                        gemo.write('vaa', vaa)
                        if date_idx == len(fd[platform])-1: del vaa
                        gemo.write('sza', sza)
                        del sza
                        gemo.write('saa', saa)
                        del saa
                        gemo.write('raa', raa)
                        del raa
                    ## end create new image

                    ## reconstruct image
                    band_data = np.vstack(band_data)
                    mask = np.vstack(mask_data)

                    ## get radiance conversion
                    slope = header[5]['slope_for_count_radiance_conversion']
                    intercept = header[5]['intercept_for_count_radiance_conversion']
                    #slope = header[5]['calibrated_slope_for_count_radiance_conversion']
                    #intercept = header[5]['calibrated_intercept_for_count_radiance_conversion']

                    ## convert to radiance
                    band_data = band_data.astype(np.float32)
                    band_data = (band_data * slope) + intercept

                    ## mask data
                    band_data[mask] = np.nan
                    band_data[sza_mask] = np.nan

                    if setu['output_lt']:
                        ds = 'Lt_{}'.format(rsrd['wave_name'][band])
                        gemo.write(ds, band_data, ds_att = ds_att)

                    ## convert to reflectance
                    scale = (np.pi * se_distance ** 2) / (ds_att['f0'] * mus)
                    band_data *= scale

                    ## mask negative TOA
                    band_data[band_data <= 0] = np.nan

                    ## write dataset
                    ds = 'rhot_{}'.format(rsrd['wave_name'][band])
                    gemo.write(ds, band_data, ds_att = ds_att)
                    if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, band_data.shape))
                    band_data = None
                gemo.close()
                gemo = None
                if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
