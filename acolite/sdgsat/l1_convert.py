## def l1_convert
## converts SDGSAT-1 KX10 MII bundle to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2023-01-03
## modifications: 2023-02-18 (QV) added merging of A and B scenes
##                2023-02-19 (QV) fixed merged scene nc_projection
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-17 (QV) use new gem NetCDF handling
##                2025-01-30 (QV) moved polygon limit
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming

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
        ## find meta and cal for given bundle
        metafiles, calfiles, imgfiles = ac.sdgsat.bundle_test(bundle)

        ## run through metafiles
        for mi, mf in enumerate(metafiles):
            cf = calfiles[mi]
            #im = imgfiles[mi]

            ## read meta and calibration
            meta = ac.sdgsat.metadata(mf)
            cal = ac.sdgsat.calibration(cf)

            ## identify sensor
            sensor = '{}_{}'.format(meta['SatelliteID'], meta['SensorID'])
            ## read rsr
            rsrd = ac.shared.rsr_dict(sensor)[sensor]

            ## get sensor specific defaults
            setd = ac.acolite.settings.parse(sensor)
            ## set sensor default if user has not specified the setting
            for k in setd:
                if k not in ac.settings['user']: setu[k] = setd[k]
            ## end set sensor specific defaults

            verbosity = setu['verbosity']
            if output is None: output = setu['output']
            if output is None: output = os.path.dirname(mf)

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
            saa = float(meta['SolarAzimuth'])
            sza = float(meta['SolarZenith'])
            vza = 5.0 ## suggested by lwk1542
            vaa = 0.0 ## azi not so important for low vza
            raa = np.abs(saa-vaa)
            while raa > 180: raa = abs(raa-360)

            dtime = dateutil.parser.parse(meta['CenterTime-Acamera'])
            isodate = dtime.isoformat()
            doy = dtime.strftime('%j')
            se_distance = ac.shared.distance_se(doy)

            ## get F0 for radiance -> reflectance computation
            f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
            f0_b = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data']*10, rsrd['rsr'])

            ## make global attributes for L1R NetCDF
            gatts = {'sensor':sensor, 'isodate':isodate, #'global_dims':global_dims,
                     'sza':sza, 'vza':vza, 'raa':raa, 'vaa': vaa, 'saa': saa,
                     'doy': doy, 'se_distance': se_distance,
                     'mus': np.cos(sza*(np.pi/180.)), 'acolite_file_type': 'L1R'}

            ## add band info to gatts
            for b in rsrd['rsr_bands']:
                gatts['{}_wave'.format(b)] = rsrd['wave_nm'][b]
                gatts['{}_name'.format(b)] = rsrd['wave_name'][b]
                gatts['{}_f0'.format(b)] = f0_b[b]

            ## output name
            oname = '{}_{}'.format(gatts['sensor'], dtime.strftime('%Y_%m_%d_%H_%M_%S'))
            if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
            ofile = '{}/{}_L1R.nc'.format(output, oname)
            gatts['oname'] = oname
            gatts['ofile'] = ofile

            ## read image projection
            dct_prj = None
            for imi, im in enumerate(imgfiles[mi]):
                print(imi, im)
                dct = ac.shared.projection_read(im)
                if (imi == 0) | (dct_prj is None):
                    sub = None
                    warp_to = None

                    ## check crop
                    if (sub is None) & (setu['limit'] is not None):
                        dct_sub = ac.shared.projection_sub(dct, setu['limit'], four_corners=True)
                        if dct_sub['out_lon']:
                            if verbosity > 1: print('Longitude limits outside {}'.format(bundle))
                            continue
                        if dct_sub['out_lat']:
                            if verbosity > 1: print('Latitude limits outside {}'.format(bundle))
                            continue
                        sub = dct_sub['sub']

                    if sub is None:
                        dct_prj = {k:dct[k] for k in dct}
                    else:
                        gatts['sub'] = sub
                        gatts['limit'] = setu['limit']
                        ## get the target NetCDF dimensions and dataset offset
                        if (warp_to is None):
                            if (setu['extend_region']): ## include part of the roi not covered by the scene
                                dct_prj = {k:dct_sub['region'][k] for k in dct_sub['region']}
                            else: ## just include roi that is covered by the scene
                                dct_prj = {k:dct_sub[k] for k in dct_sub}
                        ## end cropped
                else:
                    ## adapt full image extent
                    if (imi > 0) & (sub is None):
                        dct_prj['xrange'] = np.min((dct_prj['xrange'][0], dct['xrange'][0])), \
                                            np.max((dct_prj['xrange'][1], dct['xrange'][1]))
                        dct_prj['yrange'] = np.max((dct_prj['yrange'][0], dct['yrange'][0])), \
                                            np.min((dct_prj['yrange'][1], dct['yrange'][1]))
                        dct_prj['dimensions'] = np.round((dct_prj['yrange'][1]-dct_prj['yrange'][0]) / dct_prj['pixel_size'][1],0).astype(int), \
                                                np.round((dct_prj['xrange'][1]-dct_prj['xrange'][0]) / dct_prj['pixel_size'][0],0).astype(int)
                        dct_prj['ydim'] = dct_prj['dimensions'][0]
                        dct_prj['xdim'] = dct_prj['dimensions'][1]
            ## end run through images for projection

            ## if image projection cannot be determined or ROI is outside scene
            if dct_prj is None: continue

            ## update gatts
            gatts['scene_xrange'] = dct_prj['xrange']
            gatts['scene_yrange'] = dct_prj['yrange']
            gatts['scene_proj4_string'] = dct_prj['proj4_string']
            gatts['scene_pixel_size'] = dct_prj['pixel_size']
            gatts['scene_dims'] = dct_prj['dimensions']
            if 'zone' in dct_prj: gatts['scene_zone'] = dct_prj['zone']

            ## get projection info for netcdf
            if setu['netcdf_projection']:
                nc_projection = ac.shared.projection_netcdf(dct_prj, add_half_pixel=True)
            else:
                nc_projection = None

            ## save projection keys in gatts
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
            if (setu['output_geolocation']):
                if verbosity > 1: print('Writing geolocation lon/lat')
                lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel=True)
                gemo.write('lon', lon)
                if verbosity > 1: print('Wrote lon ({})'.format(lon.shape))
                lon = None
                gemo.write('lat', lat)
                if verbosity > 1: print('Wrote lat ({})'.format(lat.shape))
                lat = None

            ## write x/y
            if (setu['output_xy']):
                if verbosity > 1: print('Writing geolocation x/y')
                x, y = ac.shared.projection_geo(dct_prj, xy=True, add_half_pixel=True)
                gemo.write('xm', x)
                if verbosity > 1: print('Wrote xm ({})'.format(x.shape))
                x = None
                gemo.write('ym', y)
                if verbosity > 1: print('Wrote ym ({})'.format(y.shape))
                y = None

            ## run through images to load data
            for imi, im in enumerate(imgfiles[mi]):
                ## read bands
                for i, b in enumerate(cal):
                    bi = int(b)
                    band = rsrd['rsr_bands'][i]
                    print('Reading band {} from {}'.format(band, im))

                    ## read data and convert to Lt
                    md, data = ac.shared.read_band(im, idx=bi, warp_to=warp_to, gdal_meta=True)
                    nodata = data == np.uint16(0)
                    data = data.astype(float) * cal[b]['gain'] + cal[b]['bias']
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
                    f0 = gatts['{}_f0'.format(band)]/10
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
        if os.path.exists(ofile) & (new is False):

            if verbosity > 1:
                print('Conversion took {:.1f} seconds'.format(time.time()-t0))
                print('Created {}'.format(ofile))

            if setu['limit'] is not None: sub = None
            if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
