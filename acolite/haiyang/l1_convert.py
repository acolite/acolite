## def l1_convert
## converts Haiyang CZI data to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2023-02-19
## modifications: 2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2023-12-04 (QV) removed tiff support and added HDF support
##                2024-04-17 (QV) use new gem NetCDF handling
##                2025-01-30 (QV) moved polygon limit and limit buffer extension
##                2025-02-04 (QV) removed percentiles, improved settings handling
##                2025-02-10 (QV) cleaned up settings use

def l1_convert(inputfile, output = None, settings = None):

    import os, glob, dateutil.parser, datetime, time
    import numpy as np
    import h5py, scipy.ndimage
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
    ofile = None
    ofiles = []

    for fi, bundle in enumerate(inputfile):
        t0 = time.time()

        ret = ac.haiyang.bundle_test(bundle)
        image = ret['image']
        if type(image) is list: image = image[0]

        ## determine type
        bn = os.path.basename(image)
        ext = os.path.splitext(bn)[1]
        if ext == '.h5':
            image_type = 'hdf'
        elif ext == '.tiff':
            image_type = 'tiff'
        if image_type == 'tiff':
            print('Image type tiff not supported for processing {}'.format(bundle))
            print(image)
            continue

        ## parse metadata
        if verbosity > 1: print('Importing metadata from {}'.format(bundle))
        meta = ac.haiyang.metadata(ret['metadata'])
        band_names = ['L_460', 'L_560', 'L_650', 'L_825']

        ## sensor info
        sensor = '{}_{}'.format(meta['SatelliteID'], meta['SensorID']).replace('-', '')
        rsrd = ac.shared.rsr_dict(sensor)[sensor]

        ## time info
        dtime = dateutil.parser.parse(meta['CentreTime'])
        isodate = dtime.isoformat()
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)

        clat = meta['CentreLocation_Latitude']
        clon = meta['CentreLocation_Longitude']
        spos = ac.shared.sun_position(dtime, clon, clat)

        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults

        verbosity = setu['verbosity']
        sub = setu['sub']
        if output is None: output = setu['output']
        if output is None: output = os.path.dirname(bundle)

        ## get gains from settings
        gains = None
        if setu['gains']:
            if (len(setu['gains_toa']) == len(rsr_bands)) &\
               (len(setu['offsets_toa']) == len(rsr_bands)):
               gains = {}
               for bi, band in enumerate(rsr_bands):
                   gains[band] = {'gain': float(setu['gains_toa'][bi]),
                                'offset': float(setu['offsets_toa'][bi])}
            else:
                print('Use of gains requested, but provided number of gain ({}) or offset ({}) values does not match number of bands in RSR ({})'.format(len(setu['gains_toa']), len(setu['offsets_toa']), len(rsr_bands)))
                print('Provide gains in band order: {}'.format(','.join(rsr_bands)))

        ## get F0 for radiance -> reflectance computation
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        f0_b = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data']*10, rsrd['rsr'])

        warp_to = None
        if image_type == 'tiff':
            ## read image projection
            dct = ac.shared.projection_read(image)

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
            ## end cropped

        if image_type == 'hdf':
            ## open file
            f = h5py.File(image, mode='r')

            ## read atts
            h5_gatts = {a: f.attrs[a] for a in f.attrs.keys()}

            ## read lat lon
            print('Reading lat/lon')
            lat = f['Navigation Data']['Latitude'][:]/10000.
            lon = f['Navigation Data']['Longitude'][:]/10000.
            ## make subset
            if (sub is None) & (setu['limit'] is not None):
                sub = ac.shared.geolocation_sub(lat, lon, setu['limit'])
                if sub is None:
                    print('Limit {} outside of image {}'.format(setu['limit'], image))
                    continue

            if sub is not None:
                print('Using sub', sub)
                lat = lat[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
                lon = lon[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]

            ## read geometry
            print('Reading geometry')
            geometry_zoom = [1, 16]
            geometry_order = 0
            geometry_offset = 4
            vaa = scipy.ndimage.zoom(f['Navigation Data']['Satellite Azimuth Angle'][:], geometry_zoom, order =  geometry_order)[:, geometry_offset:-geometry_offset]
            if sub is not None: vaa = vaa[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
            saa = scipy.ndimage.zoom(f['Navigation Data']['Sun Azimuth Angle'][:], geometry_zoom, order = geometry_order)[:, geometry_offset:-geometry_offset]
            if sub is not None: saa = saa[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
            vza = scipy.ndimage.zoom(f['Navigation Data']['Satellite Zenith Angle'][:], geometry_zoom, order =  geometry_order)[:, geometry_offset:-geometry_offset]
            if sub is not None: vza = vza[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]
            sza = scipy.ndimage.zoom(f['Navigation Data']['Sun Zenith Angle'][:], geometry_zoom, order =  geometry_order)[:, geometry_offset:-geometry_offset]
            if sub is not None: sza = sza[sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]

            ## relative azimuth
            raa = np.abs(saa - vaa)
            raa[raa>= 180.] = np.abs(raa[raa>= 180.]-360)

        ## set global attributes
        gatts = {'sensor': sensor, 'satellite': sensor,'isodate': isodate,
                'se_distance': se_distance, 'acolite_file_type': 'L1R'}

        ## get observation geometry
        gatts['vaa'] = np.nanmean(vaa)
        gatts['saa'] = np.nanmean(saa)
        gatts['vza'] = np.nanmean(vza)
        gatts['sza'] = np.nanmean(sza)
        gatts['raa'] = np.abs(gatts['saa'] - gatts['vaa'])
        while gatts['raa'] >= 180.: gatts['raa'] = np.abs(gatts['raa']-360)
        gatts['mus'] = np.cos(gatts['sza']*(np.pi/180.))

        ## add band info to gatts
        for b in rsrd['rsr_bands']:
            gatts['{}_wave'.format(b)] = rsrd['wave_nm'][b]
            gatts['{}_name'.format(b)] = rsrd['wave_name'][b]
            gatts['{}_f0'.format(b)] = f0_b[b]

        stime = dateutil.parser.parse(gatts['isodate'])

        ## set up oname (without directory or file type) and ofile (with directory and file type)
        oname = '{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        nc_projection = None
        if image_type == 'tiff':
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
            xyr = [min(dct_prj['xrange']), min(dct_prj['yrange']),
                   max(dct_prj['xrange']), max(dct_prj['yrange']),
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
        ## end tiff projection

        ## output file
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gatts[k] for k in gatts}
        gemo.nc_projection = nc_projection

        ## write lat/lon
        if (setu['output_geolocation']):
            if verbosity > 1: print('Writing geolocation lon/lat')
            if image_type == 'tiff': lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel=True)
            gemo.write('lon', lon)
            if verbosity > 1: print('Wrote lon ({})'.format(lon.shape))
            lon = None
            gemo.write('lat', lat)
            if verbosity > 1: print('Wrote lat ({})'.format(lat.shape))
            lat = None

        ## per pixel cosine sun zenith
        mus = np.cos(np.radians(sza))

        ## write geometry
        if (setu['output_geometry']):
            if verbosity > 1: print('Writing geometry')
            ## azimuth
            gemo.write('vaa', vaa)
            if verbosity > 1: print('Wrote vaa ({})'.format(vaa.shape))
            vaa = None
            gemo.write('saa', saa)
            if verbosity > 1: print('Wrote saa ({})'.format(saa.shape))
            saa = None
            gemo.write('raa', raa)
            if verbosity > 1: print('Wrote raa ({})'.format(raa.shape))
            raa = None

            ## zenith
            gemo.write('vza', vza)
            if verbosity > 1: print('Wrote vaa ({})'.format(vza.shape))
            vza = None
            gemo.write('sza', sza)
            if verbosity > 1: print('Wrote saa ({})'.format(sza.shape))
            sza = None

        ## write x/y
        if (setu['output_xy']) & (image_type == 'tiff'):
            if verbosity > 1: print('Writing geolocation x/y')
            x, y = ac.shared.projection_geo(dct_prj, xy=True, add_half_pixel=True)
            gemo.write('xm', x)
            if verbosity > 1: print('Wrote xm ({})'.format(x.shape))
            x = None
            gemo.write('ym', y)
            if verbosity > 1: print('Wrote ym ({})'.format(y.shape))
            y = None

        ## run through bands
        for bi, band in enumerate(rsrd['rsr_bands']):
            print('Reading band {} from {}'.format(band, image))

            ## read data and convert to Lt
            if image_type == 'tiff':
                md, data = ac.shared.read_band(image, idx=bi+1, warp_to=warp_to, gdal_meta=True)
            else:
                if sub is None:
                    data = f['Geophysical Data'][band_names[bi]][:]
                else:
                    data = f['Geophysical Data'][band_names[bi]][sub[1]:sub[1]+sub[3], sub[0]:sub[0]+sub[2]]

            ## data mask
            nodata = data == 0
            print('Read band {} ({})'.format(band, data.shape))

            ds_att = {'wavelength':rsrd['wave_nm'][band]}
            if gains != None:
                ds_att['gain'] = gains[band]['gain']
                ds_att['offset'] = gains[band]['offset']
                ds_att['gains_parameter'] = setu['gains_parameter']

            if setu['output_lt']:
                ds = 'Lt_{}'.format(rsrd['wave_name'][band])
                ## write toa radiance
                gemo.write(ds, data, ds_att = ds_att)
                if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data.shape))

            ## convert to rhot
            ds = 'rhot_{}'.format(rsrd['wave_name'][band])
            data *= (np.pi * gatts['se_distance']**2) / (gatts['{}_f0'.format(band)]/100 * mus)

            ## apply gains
            if (gains != None) & (setu['gains_parameter'] == 'reflectance'):
                print('Applying gain {} and offset {} to TOA reflectance for band {}'.format(gains[band]['gain'], gains[band]['offset'], band))
                data = gains[band]['gain'] * data + gains[band]['offset']

            data[nodata] = np.nan
            if (setu['polygon_clip']): data[clip_mask] = np.nan

            ## write to netcdf file
            gemo.write(ds, data, ds_att = ds_att, replace_nan = True)
            new = False
            if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data.shape))

        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))
        gemo.close()

        if setu['limit'] is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)
        f = None

    return(ofiles, setu)
