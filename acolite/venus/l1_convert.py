## def l1_convert
## converts VENUS data to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-04-08
## modifications:  2021-12-31 (QV) new handling of settings
##                2022-01-04 (QV) added netcdf compression

def l1_convert(inputfile, output = None, settings = {},
                percentiles_compute = True,
                percentiles = (0,1,5,10,25,50,75,90,95,99,100),
                verbosity = 0):


    import os
    import dateutil.parser, time
    import numpy as np
    import acolite as ac
    import scipy.ndimage

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
        sub = None

        t0 = time.time()
        meta = ac.venus.metadata_parse(bundle)
        if meta['image_type'] != 'Reflectance':
            print('VENUS image type {} not configured'.format(meta['image_type']))
            continue

        dtime = dateutil.parser.parse(meta['isodate'])
        doy = dtime.strftime('%j')
        se_distance = ac.shared.distance_se(doy)
        isodate = dtime.isoformat()

        ## parse sensor settings
        setu = ac.acolite.settings.parse(meta['sensor'], settings=settings)
        verbosity = setu['verbosity']

        ## get other settings
        limit = setu['limit']
        output_geolocation = setu['output_geolocation']
        output_xy = setu['output_xy']
        netcdf_projection = setu['netcdf_projection']

        vname = setu['region_name']
        gains = setu['gains']
        gains_toa = setu['gains_toa']
        if output is None: output = setu['output']

        ## check if ROI polygon is given
        if setu['polylakes']:
            poly = ac.shared.polylakes(setu['polylakes_database'])
            setu['polygon_limit'] = False
        else:
            poly = setu['polygon']
        clip, clip_mask = False, None
        if poly is not None:
            if os.path.exists(poly):
                try:
                    limit = ac.shared.polygon_limit(poly)
                    if setu['polygon_limit']:
                        print('Using limit from polygon envelope: {}'.format(limit))
                    else:
                        limit = setu['limit']
                    clip = True
                except:
                    print('Failed to import polygon {}'.format(poly))

        ## read rsr
        rsrf = ac.path+'/data/RSR/{}.txt'.format(meta['sensor'])
        rsr, rsr_bands = ac.shared.rsr_read(rsrf)
        waves = np.arange(250, 2500)/1000
        waves_mu = ac.shared.rsr_convolute_dict(waves, waves, rsr)
        waves_names = {'{}'.format(b):'{:.0f}'.format(waves_mu[b]*1000) for b in waves_mu}

        ## get F0 - not stricty necessary if using reflectance
        f0 = ac.shared.f0_get()
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data'])*10, rsr)


        gatts = {'sensor':meta['sensor'],
                     'isodate':isodate, #'global_dims':global_dims,
                     'sza':meta['sza'], 'vza':meta['vza'], 'raa':meta['raa'],
                     'se_distance': se_distance,
                     'mus': np.cos(meta['sza']*(np.pi/180.)), 'acolite_file_type': 'L1R'}


        stime = dateutil.parser.parse(gatts['isodate'])
        oname = '{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'))
        if vname != '': oname+='_{}'.format(vname)
        ofile = '{}/{}_L1R.nc'.format(output, oname)

        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## add band info to gatts
        for b in rsr_bands:
            gatts['{}_wave'.format(b)] = waves_mu[b]*1000
            gatts['{}_name'.format(b)] = waves_names[b]
            gatts['{}_f0'.format(b)] = f0_b[b]

        #gatts['scene_xrange'] = meta['xrange']
        #gatts['scene_yrange'] = meta['yrange']
        #gatts['scene_pixel_size'] = meta['pixel_size']
        #gatts['scene_dims'] = meta['scene_dims']
        #gatts['scene_proj4_string'] = meta['proj4_string']
        #if 'zone' in meta: gatts['scene_zone'] = meta['zone']

        ## get projection and subset
        btag = 'B{}'.format(1)
        image_file = meta['bands'][btag]['path']
        dct = ac.shared.projection_read(image_file)
        if limit is not None:
            dct_sub = ac.shared.projection_sub(dct, limit, four_corners=True)

        ## check crop
        if (sub is None) & (limit is not None):
            dct_sub = ac.shared.projection_sub(dct, limit, four_corners=True)
            if dct_sub['out_lon']:
                if verbosity > 1: print('Longitude limits outside {}'.format(bundle))
                #continue
            if dct_sub['out_lat']:
                if verbosity > 1: print('Latitude limits outside {}'.format(bundle))
                #continue
            sub = dct_sub['sub']

        warp_to = None
        if sub is None:
            dct_prj = {k:dct[k] for k in dct}
        else:
            gatts['sub'] = sub
            gatts['limit'] = limit
            ## get the target NetCDF dimensions and dataset offset
            if (warp_to is None):
                if (extend_region): ## include part of the roi not covered by the scene
                    dct_prj = {k:dct_sub['region'][k] for k in dct_sub['region']}
                else: ## just include roi that is covered by the scene
                    dct_prj = {k:dct_sub[k] for k in dct_sub}
            ## end cropped

        pkeys = ['xrange', 'yrange', 'proj4_string', 'pixel_size', 'zone']
        for k in pkeys:
            if k in dct_prj: gatts[k] = dct_prj[k]


        ## warp settings for read_band
        ## remove extra pixel here
        xyr = [min(dct_prj['xrange']),
                min(dct_prj['yrange']),
                max(dct_prj['xrange']),
                max(dct_prj['yrange']),
                dct_prj['proj4_string']]

        res_method = 'average'
        warp_to = (dct_prj['proj4_string'], xyr, dct_prj['pixel_size'][0],dct_prj['pixel_size'][1], res_method)

        ## store scene and output dimensions
        gatts['scene_dims'] = dct['ydim'], dct['xdim']
        gatts['global_dims'] = dct_prj['dimensions']

        ## if we are clipping to a given polygon get the clip_mask here
        if clip:
            clip_mask = ac.shared.polygon_crop(dct_prj, poly, return_sub=False)
            clip_mask = clip_mask.astype(bool) == False

        ## write lat/lon
        if (output_geolocation):
            if (os.path.exists(ofile) & (not new)):
                datasets = ac.shared.nc_datasets(ofile)
            else:
                datasets = []
            if ('lat' not in datasets) or ('lon' not in datasets):
                if verbosity > 1: print('Writing geolocation lon/lat')
                lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel=True)
                ac.output.nc_write(ofile, 'lon', lon, attributes=gatts, new=new, double=True,
                                    netcdf_compression=setu['netcdf_compression'],
                                    netcdf_compression_level=setu['netcdf_compression_level'])
                lon = None
                if verbosity > 1: print('Wrote lon')
                ac.output.nc_write(ofile, 'lat', lat, double=True,
                                    netcdf_compression=setu['netcdf_compression'],
                                    netcdf_compression_level=setu['netcdf_compression_level'])
                lat = None
                if verbosity > 1: print('Wrote lat')
                new=False

        ## write x/y
        if (output_xy):
            if os.path.exists(ofile) & (not new):
                datasets = ac.shared.nc_datasets(ofile)
            else:
                datasets = []
            if ('x' not in datasets) or ('y' not in datasets):
                if verbosity > 1: print('Writing geolocation x/y')
                x, y = ac.shared.projection_geo(dct_prj, xy=True, add_half_pixel=True,
                                    netcdf_compression=setu['netcdf_compression'],
                                    netcdf_compression_level=setu['netcdf_compression_level'])
                ac.output.nc_write(ofile, 'x', x, new=new)
                x = None
                if verbosity > 1: print('Wrote x')
                ac.output.nc_write(ofile, 'y', y,
                                    netcdf_compression=setu['netcdf_compression'],
                                    netcdf_compression_level=setu['netcdf_compression_level'])
                y = None
                if verbosity > 1: print('Wrote y')
                new=False

        ## add cloud mask
        if True:
            im = None
            for k in meta['data']:
                if im is not None: continue
                if k['name'] == 'Cloud_Altitude_Grid':
                    im = k['path']
            cla = ac.shared.read_band(im, idx=1, warp_to=warp_to)
            ac.output.nc_write(ofile, 'cla', cla, new=new,
                                netcdf_compression=setu['netcdf_compression'],
                                netcdf_compression_level=setu['netcdf_compression_level'])
            cla = None
            if verbosity > 1: print('Wrote cla')

        ## read resolved geometry (only azimuths are provided?)
        if False:
            ## get solar azimuth grid
            im = None
            for k in meta['data']:
                if im is not None: continue
                if k['name'] == 'Solar_Angles_Grid':
                    im = k['path']
            ## 3000m
            tmp1 = ac.shared.read_band(im, idx=1, warp_to=warp_to)
            tmp2 = ac.shared.read_band(im, idx=2, warp_to=warp_to)
            ## 8000m
            #tmp1 = ac.shared.read_band(im, idx=3)
            #tmp2 = ac.shared.read_band(im, idx=4)
            saa = np.arctan2(tmp1,tmp2)
            saa *= 180/np.pi

            ## get view azimuth grid
            im = None
            for k in meta['data']:
                if im is not None: continue
                if k['name'] == 'Viewing_Angles_Grid':
                    im = k['path']

            ## different detectors
            for det in (1,2,3,4):
                id_x = 1 + (det-1) * 2
                id_y = 2 + (det-1) * 2
                tmp1 = ac.shared.read_band(im, idx=id_x, warp_to=warp_to)
                tmp2 = ac.shared.read_band(im, idx=id_y, warp_to=warp_to)

                vaa = np.arctan2(tmp1,tmp2)
                vaa *= 180/np.pi

        ## convert bands
        for ib, b in enumerate(rsr_bands):
            btag = 'B{}'.format(b)
            image_file = meta['bands'][btag]['path']

            ## read data
            data = ac.shared.read_band(image_file, idx=1, warp_to=warp_to)
            nodata = data == np.uint16(0)


            data = data.astype(float) / meta['REFLECTANCE_QUANTIFICATION_VALUE']
            data[nodata] = np.nan
            print(data.shape)
            ## clip to poly
            if clip: data[clip_mask] = np.nan

            ds = 'rhot_{}'.format(waves_names[b])
            ds_att = {'wavelength':waves_mu[b]*1000}
            if percentiles_compute:
                ds_att['percentiles'] = percentiles
                ds_att['percentiles_data'] = np.nanpercentile(data, percentiles)

            ## write to netcdf file
            ac.output.nc_write(ofile, ds, data, replace_nan=True, attributes=gatts, new=new, dataset_attributes = ds_att,
                            netcdf_compression=setu['netcdf_compression'],
                            netcdf_compression_level=setu['netcdf_compression_level'],
                            netcdf_compression_least_significant_digit=setu['netcdf_compression_least_significant_digit'])
            new = False
            if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data.shape))

        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))

        if limit is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
