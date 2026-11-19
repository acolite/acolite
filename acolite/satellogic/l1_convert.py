## def l1_convert
## converts Satellogic bundle to l1r NetCDF for acolite
## note: works with L1C sample data from s3://satellogic-sample-data
## written by Quinten Vanhellemont, RBINS
## 2026-11-17
## modifications: 2026-11-19 (QV) added support for L1D files

def l1_convert(inputfile, output = None, settings = None):
    import os, json, glob
    import numpy as np
    import dateutil.parser
    import acolite as ac

    from osgeo import gdal
    gdal.UseExceptions()

    ## get run/user/sensor settings
    setu = ac.acolite.settings.merge(sensor = None, settings = settings)

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if setu['verbosity'] > 1: print('Starting conversion of {} scene{}'.format(nscenes, '' if nscenes==1 else 's'))

    ## list to store output files
    ofiles = []

    ## run through inputfiles
    for bi, bundle in enumerate(inputfile):

        ## get image and metadata file
        items = ac.satellogic.bundle_test(bundle)
        if len(items) == 0: continue

        acolite_file_type = 'L1R'
        ds_base = 'rhot'
        for ii, it in enumerate(items):
            image_key = None
            if 'analytic' in items[it]: # L1C
                image_key = 'analytic'
            elif 'TOA_0' in items[it]: # L1D
                image_key = 'TOA_0'
            if image_key is None: continue

            with open(items[it]['meta'], 'r') as f:
                meta = json.load(f)

            date = meta['properties']['datetime']
            dtime = dateutil.parser.parse(date)

            if 'satl:sat_id' in meta['properties']: ## L1D
                sensor = 'NewSat'+meta['properties']['satl:sat_id'][6:]
            else: ## L1C
                sensor = 'NewSat'+meta['properties']['platform'][6:]

            vaa = meta['properties']['view:azimuth']
            #meta['properties']['view:off_nadir']
            vza = meta['properties']['view:incidence_angle']

            saa = meta['properties']['view:sun_azimuth']
            sza = 90 - meta['properties']['view:sun_elevation']
            raa = np.abs(saa - vaa)
            if raa > 180: raa = np.abs(360-raa)

            ## update settings
            setu = ac.acolite.settings.merge(sensor = sensor, settings = settings)

            ## set output directory
            if output is None: output = setu['output']
            if output is None: output = os.path.dirname(file)

            sensor_lut = '{}'.format(sensor)
            if 'sensor_lut' in setu:
                if setu['sensor_lut'] is not None:
                    sensor_lut = setu['sensor_lut']

            isodate = dtime.isoformat()

            ## set up gatts
            gatts = {}
            gatts['sensor'] = sensor
            gatts['sensor_lut'] = sensor_lut
            gatts['isodate'] = isodate
            gatts['sza'] = sza
            gatts['saa'] = saa
            gatts['vza'] = vza
            gatts['vaa'] = vaa
            gatts['raa'] = raa

            ## read RSR
            rsrd = ac.shared.rsr_dict(sensor = sensor_lut)[sensor_lut]

            ## read F0 data
            f0 = ac.shared.f0_get(f0_dataset = setu['solar_irradiance_reference'])
            f0_b = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], rsrd['rsr'])

            ## parse date
            time = dateutil.parser.parse(gatts['isodate'])
            doy = time.strftime('%j')
            se_distance = ac.shared.distance_se(doy)

            ## add to gatts
            gatts['doy'] = doy
            gatts['se_distance'] = se_distance
            gatts['acolite_file_type'] =  acolite_file_type

            ## output file name
            oname = '{}_{}'.format(gatts['sensor'], time.strftime('%Y_%m_%d_%H_%M_%S'))
            if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
            ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
            gatts['oname'] = oname
            gatts['ofile'] = ofile

            ## add band info to gatts
            for b in rsrd['rsr_bands']:
                gatts['{}_wave'.format(b)] = rsrd['wave_nm'][b]
                gatts['{}_name'.format(b)] = rsrd['wave_name'][b]
                gatts['{}_f0'.format(b)] = f0_b[b]

            ## image file path
            image_file = items[it][image_key]

            ## set up projection
            warp_to, dct_prj, sub = None, None, None
            nc_projection = None
            try:
                ## get projection from image
                dct = ac.shared.projection_read(image_file)
            except:
                if setu['verbosity'] > 1: print('Could not determine image projection')
                dct = None

            ## find crop
            if (setu['limit'] is not None) and (dct is not None):
                dct_sub = ac.shared.projection_sub(dct, setu['limit'])
                if dct_sub['out_lon']:
                    if setu['verbosity'] > 1: print('Longitude limits outside {}'.format(bundle))
                    continue
                if dct_sub['out_lat']:
                    if setu['verbosity'] > 1: print('Latitude limits outside {}'.format(bundle))
                    continue
                sub = dct_sub['sub']

            if dct is not None:
                if sub is None:
                    dct_prj = {k:dct[k] for k in dct}
                else:
                    dct_prj = {k:dct_sub[k] for k in dct_sub}
                    xyr = [min(dct_prj['xrange']),
                           min(dct_prj['yrange']),
                           max(dct_prj['xrange']),
                           max(dct_prj['yrange']),
                           dct_prj['proj4_string']]
                    res_method = 'near'
                    warp_to = (dct_prj['proj4_string'], xyr, dct_prj['pixel_size'][0],dct_prj['pixel_size'][1], res_method)

                ## set up nc_projection
                nc_projection = ac.shared.projection_netcdf(dct_prj)
            else:
                ret = ac.shared.warp_and_merge(items[it][image_key], output = output,
                                               limit = setu['limit'], resolution = setu['default_projection_resolution'])
                if len(ret) == 3:
                    image_file = ret[0]
                    dct = ret[1]
                    dct_limit = ret[2]
                    dct_prj = dct

            ## for radiance conversion
            #mus = np.cos(gatts['sza']*(np.pi/180))
            #muv = np.cos(gatts['vza']*(np.pi/180))

            ## set up output gem
            gemo = ac.gem.gem(ofile, new = True)
            gemo.gatts = {k: gatts[k] for k in gatts}
            gemo.nc_projection = nc_projection

            ## output geolocation
            if dct_prj is not None:
                if setu['verbosity'] > 1: print('Computing and writing lat/lon')
                ## offset half pixels to compute center pixel lat/lon
                dct_prj['xrange'] = dct_prj['xrange'][0]+dct_prj['pixel_size'][0]/2, dct_prj['xrange'][1]-dct_prj['pixel_size'][0]/2
                dct_prj['yrange'] = dct_prj['yrange'][0]+dct_prj['pixel_size'][1]/2, dct_prj['yrange'][1]-dct_prj['pixel_size'][1]/2
                ## compute lat/lon
                lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel = False)
                print(lat.shape)
                gemo.write('lat', lat)
                lat = None
                gemo.write('lon', lon)
                lon = None

            ## run through bands
            for bi, band in enumerate(rsrd['rsr_bands']):

                ds_att = {'wavelength': rsrd['wave_nm'][band], 'wave_name': rsrd['wave_name'][band] }
                wave = ds_att['wave_name']
                ds = '{}_{}'.format(ds_base, ds_att['wave_name'])
                if setu['verbosity'] > 2: print('Computing {} for {}'.format(ds, gatts['oname']))

                cdata = ac.shared.read_band(image_file, idx = 1 + bi, sub = sub, warp_to = warp_to)
                cdata = cdata.astype(np.float32) * 0.0001

                ## write toa radiance
                if setu['output_lt']: print('Lt writing not configured for Satellogic')

                ## write reflectance
                gemo.write(ds, cdata, ds_att = ds_att)
                cdata = None

            ## update attributes
            gemo.gatts = {k: gatts[k] for k in gatts}
            gemo.gatts_update()
            gemo.close()

            ofiles.append(ofile)
    return(ofiles, setu)
