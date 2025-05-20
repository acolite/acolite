## def l1_convert
## converts DEIMOS2 bundle to l1r NetCDF
## written by Quinten Vanhellemont, RBINS
## 2024-05-16
## modifications: 2025-01-30 (QV) moved polygon limit
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming
##                2025-05-15 (QV) added Deimos1
##                2025-05-20 (QV) added rsr_version

def l1_convert(inputfile, output = None, settings = None):
    import os
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
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    new = True
    warp_to = None

    ofile = None
    ofiles = []

    for bundle in inputfile:
        t0 = time.time()

        ## get info
        images = ac.deimos.bundle_test(bundle)
        meta = ac.deimos.metadata(images['MS']['metadata'])

        if 'PAN' in images:
            meta_pan = ac.deimos.metadata(images['PAN']['metadata'])
        else:
            meta_pan = None

        if (meta['MISSION'] == 'Deimos 2') & (meta['INSTRUMENT'] == 'HiRAIS'):
            sensor = 'DEIMOS2_HiRAIS'
        elif (meta['MISSION'] == 'DEIMOS') & (meta['MISSION_INDEX'] == '1'):
            sensor = 'DEIMOS1_SLIM6'
        else:
            continue

        print('Warning: Experimental support for {}'.format(sensor))
        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults

        verbosity = setu['verbosity']
        if output is None: output = setu['output']
        if output is None: output = os.path.dirname(bundle)

        ## read rsr
        if setu['rsr_version'] is not None:
            if setu['rsr_version'] == 'square':
                print('Warning: Using square RSR for {}'.format(sensor))
            sensor_version = '{}_{}'.format(sensor, setu['rsr_version'])
        else:
            sensor_version = '{}'.format(sensor)
        rsrd = ac.shared.rsr_dict(sensor_version)[sensor_version]

        gains = None
        if setu['gains']:
            if (len(setu['gains_toa']) == len(rsrd['rsr_bands'])) &\
                (len(setu['offsets_toa']) == len(rsrd['rsr_bands'])):
                gains = {}
                for bi, band in enumerate(rsrd['rsr_bands']):
                    gains[band] = {'gain': float(setu['gains_toa'][bi]),
                                'offset': float(setu['offsets_toa'][bi])}
            else:
                print('Use of gains requested, but provided number of gain ({}) or offset ({}) values does not match number of bands in RSR ({})'.format(len(setu['gains_toa']), len(setu['offsets_toa']), len(rsr_bands)))
                print('Provide gains in band order: {}'.format(','.join(rsrd['rsr_bands'])))

        ## geometry
        saa = float(meta['SUN_AZIMUTH'])
        sza = 90 - float(meta['SUN_ELEVATION'])
        vza = abs(float(meta['VIEWING_ANGLE']))
        vaa = 0

        raa = np.abs(saa-vaa)
        while raa > 180: raa = abs(raa-360)

        if 'IMAGING_TIME' in meta:
            dtime = dateutil.parser.parse(meta['IMAGING_TIME'])
        else:
            dtime = dateutil.parser.parse(meta['START_TIME'])
        isodate = dtime.isoformat()
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)

        ## get F0 for radiance -> reflectance computation
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        f0_b = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], rsrd['rsr'])

        ## make global attributes for L1R NetCDF
        gatts = {'sensor':sensor, 'isodate':isodate,
                 'sza':sza, 'vza':vza, 'raa':raa, 'vaa': vaa, 'saa': saa,
                 'doy': doy, 'se_distance': se_distance,
                 'mus': np.cos(sza*(np.pi/180.)), 'acolite_file_type': 'L1R'}

        ## add band info to gatts
        for b in rsrd['rsr_bands']:
            gatts['{}_wave'.format(b)] = rsrd['wave_nm'][b]
            gatts['{}_name'.format(b)] = rsrd['wave_name'][b]
            gatts['{}_f0'.format(b)] = f0_b[b]

        ## output file name
        oname = '{}_{}'.format(gatts['sensor'], dtime.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## identify projection
        try:
            dct_prj = ac.shared.projection_read(images['MS']['image'])
        except:
            dct_prj = None

        ## test subset
        sub = None
        if (sub is None) & (setu['limit'] is not None):
            dct_sub = ac.shared.projection_sub(dct_prj, setu['limit'], four_corners=True)
            if dct_sub['out_lon']:
                if verbosity > 1: print('Longitude limits outside {}'.format(bundle))
                continue
            if dct_sub['out_lat']:
                if verbosity > 1: print('Latitude limits outside {}'.format(bundle))
                continue
            sub = dct_sub['sub']
            dct_prj = dct_sub

        ## add projection keys to gatts
        nc_projection = None
        if dct_prj is not None:
            ## get projection info for netcdf
            if setu['netcdf_projection']:
                nc_projection = ac.shared.projection_netcdf(dct_prj, add_half_pixel=True)

            pkeys = ['xrange', 'yrange', 'proj4_string', 'pixel_size', 'zone']
            for k in pkeys:
                if k in dct_prj: gatts[k] = dct_prj[k]

            ## warp settings for read_band
            xyr = [min(dct_prj['xrange']),
                    min(dct_prj['yrange']),
                    max(dct_prj['xrange']),
                    max(dct_prj['yrange']),
                    dct_prj['proj4_string']]

            res_method = 'average'
            warp_to = (dct_prj['proj4_string'], xyr, dct_prj['pixel_size'][0],dct_prj['pixel_size'][1], res_method)

            ## store scene and output dimensions
            gatts['scene_dims'] = dct_prj['ydim'], dct_prj['xdim']
            gatts['global_dims'] = dct_prj['dimensions']

            ## if we are clipping to a given polygon get the clip_mask here
            if setu['polygon_clip']:
                clip_mask = ac.shared.polygon_crop(dct_prj, setu['polygon'], return_sub=False)
                clip_mask = clip_mask.astype(bool) == False

        ## output file
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gatts[k] for k in gatts}
        gemo.nc_projection = nc_projection
        new = False

        ## write lat/lon
        if (setu['output_geolocation']) & (dct_prj is not None):
            if verbosity > 1: print('Writing geolocation lon/lat')
            lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel=True)
            gemo.write('lon', lon)
            if verbosity > 1: print('Wrote lon ({})'.format(lon.shape))
            lon = None
            gemo.write('lat', lat)
            if verbosity > 1: print('Wrote lat ({})'.format(lat.shape))
            lat = None

        ## write x/y
        if (setu['output_xy'])  & (dct_prj is not None):
            if verbosity > 1: print('Writing geolocation x/y')
            x, y = ac.shared.projection_geo(dct_prj, xy=True, add_half_pixel=True)
            gemo.write('xm', x)
            if verbosity > 1: print('Wrote xm ({})'.format(x.shape))
            x = None
            gemo.write('ym', y)
            if verbosity > 1: print('Wrote ym ({})'.format(y.shape))
            y = None

        ## read bands
        for i, band in enumerate(rsrd['rsr_bands']):
            if band == 'PAN':
                if 'PAN' not in images: continue
                im = images['PAN']['image']
                gain = meta_pan['bands']['PAN']['PHYSICAL_GAIN']
                bias = meta_pan['bands']['PAN']['PHYSICAL_BIAS']
                if 'ESUN' in meta_pan['bands']['PAN']:
                    esun = meta_pan['bands']['PAN']['ESUN']
                else:
                    esun = gatts['{}_f0'.format(band)]
                bi = meta_pan['bands']['PAN']['BAND_INDEX']
            else:
                im = images['MS']['image']
                gain = meta['bands'][band]['PHYSICAL_GAIN']
                bias = meta['bands'][band]['PHYSICAL_BIAS']
                if 'ESUN' in meta['bands'][band]:
                    esun = meta['bands'][band]['ESUN']
                else:
                    esun = gatts['{}_f0'.format(band)]
                bi = meta['bands'][band]['BAND_INDEX']
                print(bi, band, esun)
            print('Reading band {} from {}'.format(band, im))

            ## read data and convert to Lt
            md, data = ac.shared.read_band(im, idx=bi, warp_to=warp_to, gdal_meta=True)
            nodata = data == np.uint16(0)
            data = data.astype(float) * gain + bias
            print('Read band {} ({})'.format(band, data.shape))

            ## apply gains
            if (gains != None) & (setu['gains_parameter'] == 'radiance'):
                print('Applying gain {} and offset {} to TOA radiance for band {}'.format(gains[band]['gain'], gains[band]['offset'], band))
                data = gains[band]['gain'] * data + gains[band]['offset']

            ds_att = {'wavelength':rsrd['wave_nm'][band]}
            if gains != None:
                ds_att['gain'] = gains[band]['gain']
                ds_att['offset'] = gains[band]['offset']
                ds_att['gains_parameter'] = setu['gains_parameter']

            if setu['output_lt']:
                ds = 'Lt_{}'.format(rsrd['wave_name'][band])
                ## write toa radiance
                gemo.write(ds, data, ds_att = ds_att)

            ## convert to rhot
            ds = 'rhot_{}'.format(rsrd['wave_name'][band])
            #f0 = gatts['{}_f0'.format(band)]/10
            f0 = esun * 1.0
            data *= (np.pi * gatts['se_distance']**2) / (f0 * np.nanmean(gatts['mus']))

            ## apply gains
            if (gains != None) & (setu['gains_parameter'] == 'reflectance'):
                print('Applying gain {} and offset {} to TOA reflectance for band {}'.format(gains[band]['gain'], gains[band]['offset'], band))
                data = gains[band]['gain'] * data + gains[band]['offset']

            data[nodata] = np.nan
            if (setu['polygon_clip']): data[clip_mask] = np.nan

            ## write to netcdf file
            gemo.write(ds, data, ds_att = ds_att, replace_nan = True)
            if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data.shape))

        ## close file
        gemo.close()

        ## check if file exists and if it was created now
        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))

        if setu['limit'] is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
