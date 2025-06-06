## def l1_convert
## convert ABI image to L1R ACOLITE NetCDF file
##
## written by Quinten Vanhellemont, RBINS
## 2025-06-05
## modifications:

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
    for bundle in inputfile:
        ## find goes files
        bundle_files = ac.goes.bundle_test(bundle)

        for satellite in bundle_files:
            platform = satellite.replace('G', 'GOES')
            sensor = '{}_ABI'.format(platform)
            rsrd = ac.shared.rsr_dict(sensor)[sensor]

            bands = rsrd['rsr_bands']
            band_fname = ['C{}'.format(bi+1)  for bi, band in enumerate(bands)]

            ## get sensor specific defaults
            setd = ac.acolite.settings.parse(sensor)
            ## set sensor default if user has not specified the setting
            for k in setd:
                if k not in ac.settings['user']: setu[k] = setd[k]
            ## end set sensor specific defaults
            if output is None: output = setu['output']

            ## get F0 for radiance -> reflectance computation
            #f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
            #f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], rsrd['rsr'])

            ## target resolution
            ssd = float(setu['abi_target_res'])
            if ssd not in [0.5, 1.0, 2.0]:
                print('abi_target_res={} not supported'.format(setu['abi_target_res']))
                continue

            ## get reference band for requested resolution
            reference_band = {0.5: 'CH064', 1.0: 'CH047', 2.0: 'CH225'}[ssd]
            bands_res = [1.0, 0.5, 1.0, 2.0, 1.0, 2.0]

            for product in bundle_files[satellite]:
                #RadF Full
                #RadC CONUS
                #RadM1 Mesoscale Area1
                #RadM2 Mesoscale Area2

                ## region sizes y, x
                if product == 'RadF':
                    if ssd == 0.5:
                        shape = 21696, 21696
                    elif ssd == 1:
                        shape = 10848, 10848
                    elif ssd == 2:
                        shape = 5424, 5424
                elif product == 'RadC':
                    if ssd == 0.5:
                        shape = 6000, 10000
                    elif ssd == 1:
                        shape = 3000, 5000
                    elif ssd == 2:
                        shape = 1500, 2500
                else:
                    if ssd == 0.5:
                        shape = 2000, 2000
                    elif ssd == 1:
                        shape = 1000, 1000
                    elif ssd == 2:
                        shape = 500, 500


                for start_date in bundle_files[satellite][product]:
                    band_files = bundle_files[satellite][product][start_date]
                    band_keys = list(band_files.keys())

                    print(start_date, len(band_files))

                    ## track output gatts
                    gatts = {}

                    ## track if we have loaded current projection and geometry
                    cur_geometry = False

                    ## set up geolocation and geometry
                    if cur_geometry == False:
                        for bi, band in enumerate(bands):
                            if band != reference_band: continue
                            reference_band_key = [k for k in band_keys if k.endswith('C{}'.format(str(bi+1).zfill(2)))]
                            if len(reference_band_key) != 1:
                                print('Geometry reference band {} not in band files'.format(reference_band))
                                continue
                            else:
                                reference_band_key = reference_band_key[0]

                        sub = None

                        ## open file
                        gemi = ac.gem.gem(band_files[reference_band_key]['path'])
                        goes_imager_projection = gemi.atts('goes_imager_projection')
                        lon_0 = goes_imager_projection['longitude_of_projection_origin']
                        gemi.close()
                        gemi = None

                        ## do we need to compute geolocation and geometry
                        compute_geol, compute_geom = True, True

                        ## find geolocation and geometry files for this product
                        geol, geom = None, None
                        if product in ['RadF', 'RadC']:
                            geol = ac.config['data_dir'] + '/GEO/{}/{}_lonlat_{}_{:.1f}km.nc'.format('ABI', product, lon_0, ssd)
                            geom = ac.config['data_dir'] + '/GEO/{}/{}_vaavza_{}_{:.1f}km.nc'.format('ABI', product, lon_0, ssd)

                        ## load geolocation
                        if geol is not None:
                            if os.path.exists(geol):
                                print('Reading geolocation from {}'.format(geol))
                                lon = ac.shared.nc_data(geol, 'lon')
                                lat = ac.shared.nc_data(geol, 'lat')
                                compute_geol = False

                        ## load geometry
                        if geom is not None:
                            if os.path.exists(geom):
                                print('Reading geometry from {}'.format(geom))
                                vaa = ac.shared.nc_data(geom, 'vaa')
                                vza = ac.shared.nc_data(geom, 'vza')
                                compute_geom = False

                        ## compute geolocation
                        if compute_geol:
                            print('Computing geolocation for {}'.format(product))
                            lon, lat = ac.goes.geolocation(band_files[reference_band_key]['path'])

                            ## write geolocation
                            if geol is not None:
                                if not os.path.exists(geol):
                                    print('Saving geolocation to {}'.format(geol))
                                    gemo = ac.gem.gem(geol, new = True)
                                    gemo.write('lon', lon)
                                    gemo.write('lat', lat)
                                    gemo.close()
                                    gemo = None
                                    print('Saved geolocation to {}'.format(geol))

                        ## compute geometry
                        if compute_geom:
                            print('Computing geometry for {}'.format(product))
                            vaa, vza = ac.goes.geometry(lon, lat, lon_0 = lon_0)

                            ## write geometry
                            if geom is not None:
                                if not os.path.exists(geom):
                                    print('Saving geometry to {}'.format(geom))
                                    gemo = ac.gem.gem(geom, new = True)
                                    gemo.write('vaa', vaa)
                                    gemo.write('vza', vza)
                                    gemo.close()
                                    gemo = None
                                    print('Saved geometry to {}'.format(geom))

                        ## find image subset
                        sub = None
                        if setu['sub'] is not None:
                            sub = setu['sub']
                        elif setu['limit'] is not None:
                            sub = ac.shared.geolocation_sub(lat, lon, setu['limit'])
                            if sub is None:
                                print('Limit {} not in current product {}.'.format(setu['limit'], product))
                                #continue
                        elif product in ['RadC', 'RadF']:
                            print('Warning: running without subsetting {} image.'.format(product))
                            print('It is recommended using a subset sub=x, y, nx, ny or limit=S, W, N, E.')
                            #return

                        ## subset geolocation and geometry
                        column_range, row_range = None, None
                        if sub is not None:
                            ## offset to multiple of 4 pixels
                            sub_mod = [s % 4 for s in sub]
                            sub = [sub[0] - sub_mod[0], sub[1] - sub_mod[1], sub[2] + sub_mod[2], sub[3] + sub_mod[3]]

                            column_range = sub[0], sub[0]+sub[2]
                            row_range = sub[1], sub[1]+sub[3]
                            ## subset geolocation
                            lon = lon[row_range[0]:row_range[1], column_range[0]:column_range[1]]
                            lat = lat[row_range[0]:row_range[1], column_range[0]:column_range[1]]
                            ## subset geometry
                            vaa = vaa[row_range[0]:row_range[1], column_range[0]:column_range[1]]
                            vza = vza[row_range[0]:row_range[1], column_range[0]:column_range[1]]

                        ## track target shape
                        target_shape = lon.shape
                        cur_geometry = True

                    ## run through bands
                    for bi, band in enumerate(bands):
                        band_key = [k for k in band_keys if k.endswith('C{}'.format(str(bi+1).zfill(2)))]
                        if len(band_key) != 1:
                            print('Band {} not in band files'.format(band))
                            continue
                        else:
                            band_key = band_key[0]

                        ## set up gatts and create new file
                        if len(gatts) == 0:
                            ## get observation time
                            start_time = dateutil.parser.parse(band_files[band_key]['gatts']['time_coverage_start'])
                            end_time = dateutil.parser.parse(band_files[band_key]['gatts']['time_coverage_end'])
                            dt = start_time + (end_time - start_time)

                            ## get up gatts
                            gatts['sensor'] = sensor
                            gatts['isodate'] = dt.isoformat()
                            gatts['acolite_file_type'] = 'L1R'

                            ## output name
                            oname = '{}_{}'.format(gatts['sensor'], dt.strftime('%Y_%m_%d_%H_%M_%S'))
                            if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
                            ofile = '{}/{}_{}'.format(output, oname, gatts['acolite_file_type'])
                            #ofile_hrv = '{}/{}_{}_HRV.nc'.format(output, oname, gatts['acolite_file_type'])
                            ofile += '_{:.1f}km.nc'.format(ssd)
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
                        ## end create new image

                        factor = 1
                        if bands_res[bi] != ssd:
                            factor = bands_res[bi]/ssd

                        band_sub = None
                        if sub is not None:
                            if factor == 1.0:
                                band_sub = [s for s in sub]
                            else:
                                band_sub = [int(np.round(s / factor)) for s in sub]

                        ## output dataset attributes
                        ds_att = {'wavelength':rsrd['wave_nm'][band]}#, 'f0': f0d[band]}

                        print(band, band_key, ds_att)
                        print(band_sub)
                        gemi = ac.gem.gem(band_files[band_key]['path'])
                        ## get some band info
                        ds_att['band_wavelength'] = gemi.data('band_wavelength').flatten()[0]
                        ds_att['esun']  = gemi.data('esun').flatten()[0]
                        ## read data
                        d, a = gemi.data('Rad', sub = band_sub, attributes = True)
                        gemi.close()
                        gemi = None

                        print(d.shape, target_shape, factor)
                        ## resample
                        if factor != 1.0:
                            #factor = d.shape[0] / target_shape[0]
                            #factor1 = d.shape[1] / target_shape[1]
                            d = scipy.ndimage.zoom(d, zoom=factor, order=1)

                        if setu['output_lt']:
                            ds = 'Lt_{}'.format(rsrd['wave_name'][band])
                            gemo.write(ds, d, ds_att = a)

                        ## convert to reflectance
                        scale = (np.pi * se_distance ** 2) / (ds_att['esun'] * mus)
                        d *= scale

                        ## mask negative TOA
                        d[d <= 0] = np.nan

                        ## write dataset
                        ds = 'rhot_{}'.format(rsrd['wave_name'][band])
                        gemo.write(ds, d, ds_att = ds_att)
                        if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, d.shape))
                        d = None
                    gemo.close()
                    gemo = None
                if ofile not in ofiles: ofiles.append(ofile)
    return(ofiles, setu)
