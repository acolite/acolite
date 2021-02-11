## def l1_convert
## converts landsat bundle data to l1r NetCDF for acolite-gen
## written by Quinten Vanhellemont, RBINS
## 2021-02-05
## modifications: 2021-02-07 (QV) added support for tile merging and extending the size of the output file to the requested limit
##                2021-02-09 (QV) added cross zone warping support (test)
##                2021-02-10 (QV) fixed cross zone warping and added support for full tile warping
##                2021-02-11 (QV) added checks for merging tiles of the same sensor and close in time

def l1_convert(inputfile, output=None,
                limit=None, sub=None,
                output_pan = True,
                output_pan_ms = True,
                output_thermal = True,

                output_geometry = True,
                output_geolocation = True,
                output_xy = False,
                usgs_reflectance = True,

                percentiles_compute = True,
                percentiles = (0,1,5,10,25,50,75,90,95,99,100),

                merge_tiles = False,
                merge_zones = False,
                extend_region = False,

                check_sensor = True,
                check_time = True,
                max_merge_time = 600, # seconds

                verbosity = 0, vname = ''):

    import os, glob, dateutil, time
    import acolite as ac
    import scipy.ndimage
    import numpy as np
    t0 = time.time()

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    ## check if merging settings make sense
    if (limit is None) & (merge_tiles):
        if verbosity > 0: print("Merging tiles not supported without ROI limit")
        merge_tiles = False
    if merge_tiles:
        merge_zones = True
        extend_region = True

    new = True
    new_pan = True
    warp_to = None
    warp_to_pan = None

    ofile = None
    ofiles = []

    for bundle in inputfile:
        if output is None: output = os.path.dirname(bundle)
        if verbosity > 1: print('Starting conversion of {}'.format(bundle))

        mtl = glob.glob('{}/{}'.format(bundle, '*MTL.txt'))
        if len(mtl) == 0:
            if verbosity > 0: print('No metadata file found for {}'.format(bundle))
            continue
        else:
            mtl = mtl[0]

        ## read landsat metadata and check files
        if verbosity > 1: print('Importing metadata from {}'.format(bundle))
        meta = ac.landsat.metadata_read(mtl)
        fmeta = ac.landsat.metadata_bands(bundle, meta)

        ## get relevant data from meta
        if 'PRODUCT_CONTENTS' in meta: ## COLL2
            pk = 'IMAGE_ATTRIBUTES'
            ik = 'IMAGE_ATTRIBUTES'
            rk = 'PROJECTION_ATTRIBUTES'
        elif 'PRODUCT_METADATA' in meta: ## COLL1
            pk = 'PRODUCT_METADATA'
            ik = 'IMAGE_ATTRIBUTES'
            rk = 'PRODUCT_METADATA'
        spacecraft_id = meta[pk]['SPACECRAFT_ID']
        sensor_id = meta[pk]['SENSOR_ID']
        path = meta[pk]['WRS_PATH']
        row = meta[pk]['WRS_ROW']
        isodate = meta[pk]['DATE_ACQUIRED']+'T'+meta[pk]['SCENE_CENTER_TIME']
        global_dims = int(meta[rk]['REFLECTIVE_LINES']), int(meta[rk]['REFLECTIVE_SAMPLES'])

        ## some hard coded info
        sat = 'L{}'.format(spacecraft_id[-1])
        if sensor_id in ['TM']:# Landsat 5
            sen = 'TM'
            pan_bands = []
            thermal_bands = ['6']
        elif sensor_id in ['ETM']:# Landsat 7
            sen = 'ETM'
            pan_bands = ['8']
            thermal_bands = ['6_VCID_1', '6_VCID_2']
        elif sensor_id in ['OLI', 'OLI_TIRS']:# Landsat 8
            sen = 'OLI'
            pan_bands = ['8']
            thermal_bands = ['10', '11']
        else:
            print(spacecraft_id, sensor_id)
            print('Not configured')
            continue

        sensor = '{}_{}'.format(sat,sen)

        ## scene average geometry
        vza = 0
        sza = 90-float(meta[ik]['SUN_ELEVATION'])
        raa = float(meta[ik]['SUN_AZIMUTH'])
        se_distance = float(meta[ik]['EARTH_SUN_DISTANCE'])

        ## read rsr
        rsrf = ac.path+'/data/RSR/{}.txt'.format(sensor)
        rsr, rsr_bands = ac.shared.rsr_read(rsrf)
        waves = np.arange(250, 2500)/1000
        waves_mu = ac.shared.rsr_convolute_dict(waves, waves, rsr)
        waves_names = {'{}'.format(b):'{:.0f}'.format(waves_mu[b]*1000) for b in waves_mu}

        ## get F0 - not stricty necessary if using USGS reflectance
        f0 = ac.shared.f0_get()
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data'])*10, rsr)

        ## make global attributes for L1R NetCDF
        gatts = {'sensor':sensor, 'isodate':isodate, 'global_dims':global_dims,
                 'sza':sza, 'vza':vza, 'raa':raa, 'se_distance': se_distance,
                 'mus': np.cos(sza*(np.pi/180.)), 'wrs_path': path, 'wrs_row': row}
        if merge_tiles:
            gatts['tile_code'] = 'merged'
        else:
            gatts['tile_code'] = '{}{}'.format(gatts['wrs_path'].zfill(3),gatts['wrs_row'].zfill(3))

        stime = dateutil.parser.parse(gatts['isodate'])
        oname = '{}_{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'), gatts['tile_code'])
        if vname != '': oname+='_{}'.format(vname)

        ## output file information
        if (merge_tiles is False) | (ofile is None):
            ofile = '{}/{}_L1R.nc'.format(output, oname)
            gatts['oname'] = oname
            gatts['ofile'] = ofile
        elif (merge_tiles) & (ofile is None):
            ofile = '{}/{}_L1R.nc'.format(output, oname)
            gatts['oname'] = oname
            gatts['ofile'] = ofile

        ## check if we should merge these tiles
        if (merge_tiles) & (not new) & (os.path.exists(ofile)):
                fgatts = ac.shared.nc_gatts(ofile)
                if (check_sensor) & (fgatts['sensor'] != gatts['sensor']):
                    print('Sensors do not match, skipping {}'.format(bundle))
                    continue
                if check_time:
                    tdiff = dateutil.parser.parse(fgatts['isodate'])-dateutil.parser.parse(gatts['isodate'])
                    tdiff = abs(tdiff.days*86400 + tdiff.seconds)
                    if (tdiff > max_merge_time):
                        print('Time difference too large, skipping {}'.format(bundle))
                        continue

        print(ofile)

        ## add band info to gatts
        for b in rsr_bands:
            gatts['{}_wave'.format(b)] = waves_mu[b]*1000
            gatts['{}_name'.format(b)] = waves_names[b]
            gatts['{}_f0'.format(b)] = f0_b[b]
            if b in fmeta:
                fmeta[b]['f0'] = f0_b[b]
                fmeta[b]['se_distance'] = gatts['se_distance']

        ## get scene projection and extent
        dct = ac.landsat.projection(meta)

        ## full scene
        gatts['scene_xrange'] = dct['xrange']
        gatts['scene_yrange'] = dct['yrange']
        gatts['scene_proj4_string'] = dct['proj4_string']
        gatts['scene_pixel_size'] = dct['pixel_size']
        gatts['scene_dims'] = dct['dimensions']
        if 'zone' in dct: gatts['scene_zone'] = dct['zone']

        if (sub is None) & (limit is not None):
            dct_sub = ac.shared.projection_sub(dct, limit, four_corners=True)
            if dct_sub['out_lon']:
                if verbosity > 1: print('Longitude limits outside {}'.format(bundle))
                continue
            if dct_sub['out_lat']:
                if verbosity > 1: print('Latitude limits outside {}'.format(bundle))
                continue
            sub = dct_sub['sub']
        else:
            if extend_region:
                print("Can't extend region if no ROI limits given")
                extend_region = False

        if ((merge_tiles is False) & (merge_zones is False)): warp_to = None
        if sub is None:
            sub_pan = None
            if ((merge_zones) & (warp_to is not None)):
                if dct_prj != dct: ## target projection differs from this tile, need to set bounds
                    if dct['proj4_string'] != dct_prj['proj4_string']:
                        ## if the prj does not match, project current scene bounds to lat/lon
                        lonr, latr = dct['p'](dct['xrange'], dct['yrange'], inverse=True)
                        ## then to target projection
                        xrange_raw, yrange_raw = dct_prj['p'](lonr, (latr[1], latr[0]))
                        ## fix to nearest full pixel
                        pixel_size = dct_prj['pixel_size']
                        dct_prj['xrange'] = [xrange_raw[0] - (xrange_raw[0] % pixel_size[0]), xrange_raw[1]+pixel_size[0]-(xrange_raw[1] % pixel_size[0])]
                        dct_prj['yrange'] = [yrange_raw[1]+pixel_size[1]-(yrange_raw[1] % pixel_size[1]), yrange_raw[0] - (yrange_raw[0] % pixel_size[1])]
                        ## need to add new dimensions
                        dct_prj['xdim'] = int((dct_prj['xrange'][1]-dct_prj['xrange'][0])/pixel_size[0])+1
                        dct_prj['ydim'] = int((dct_prj['yrange'][1]-dct_prj['yrange'][0])/pixel_size[1])+1
                        dct_prj['dimensions'] = [dct_prj['xdim'], dct_prj['ydim']]
                    else:
                        ## if the projection matches just use the current scene projection
                        dct_prj = {k:dct[k] for k in dct}
            elif (warp_to is None):
                dct_prj = {k:dct[k] for k in dct}
        else:
            pan_dims = sub[3]*2, sub[2]*2
            sub_pan = [s*2 for s in sub]
            gatts['sub'] = sub
            gatts['pan_sub'] = sub_pan
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

        ## else use the projection info in dct_prj
        xyr = [min(dct_prj['xrange'])-dct_prj['pixel_size'][0]/2,min(dct_prj['yrange']),
               max(dct_prj['xrange']),max(dct_prj['yrange'])-dct_prj['pixel_size'][1]/2,
               dct_prj['proj4_string']]
        xyr_pan = [min(dct_prj['xrange'])-dct_prj['pixel_size'][0],min(dct_prj['yrange']),
                   max(dct_prj['xrange']),max(dct_prj['yrange'])-dct_prj['pixel_size'][1],
                   dct_prj['proj4_string']]

        ## warp settings for read_band
        warp_to = (dct_prj['proj4_string'], xyr, dct_prj['pixel_size'][0],dct_prj['pixel_size'][1],'near')
        warp_to_pan = (dct_prj['proj4_string'], xyr_pan, dct_prj['pixel_size'][0]/2,dct_prj['pixel_size'][1]/2,'near')

        ## store scene and output dimensions
        gatts['scene_dims'] = dct['ydim'], dct['xdim']
        gatts['global_dims'] = dct_prj['dimensions']
        gatts['pan_dims'] =  dct_prj['dimensions'][0]*2, dct_prj['dimensions'][1]*2

        ## new file for every bundle if not merging
        if (merge_tiles is False):
            new = True
            new_pan = True

        ## start the conversion
        ## write geometry
        if ('VAA' in fmeta) & ('SAA' in fmeta) & ('VZA' in fmeta) & ('SZA' in fmeta):
            if verbosity > 1: print('Reading per pixel geometry')
            sza = ac.shared.read_band(fmeta['SZA']['FILE'], sub=sub, warp_to=warp_to).astype(np.float32)/100
            mus = np.cos(sza*(np.pi/180.)) ## per pixel cos sun zenith
            if (output_geometry):
                saa = ac.shared.read_band(fmeta['SAA']['FILE'], sub=sub, warp_to=warp_to).astype(np.float32)/100
                vza = ac.shared.read_band(fmeta['VZA']['FILE'], sub=sub, warp_to=warp_to).astype(np.float32)/100
                vaa = ac.shared.read_band(fmeta['VAA']['FILE'], sub=sub, warp_to=warp_to).astype(np.float32)/100
                mask = (vaa == 0) * (vza == 0) * (saa == 0) * (sza == 0)
                vza[mask] = np.nan
                sza[mask] = np.nan
                raa = (saa-vaa)

                tmp = np.where(raa>180)
                raa[tmp]=np.abs(raa[tmp] - 360)
                raa[mask] = np.nan

                vaa = None
                saa = None
                mask = None
                ac.output.nc_write(ofile, 'raa', raa, replace_nan=True, attributes=gatts, new=new)
                if verbosity > 1: print('Wrote raa')
                new = False
                ac.output.nc_write(ofile, 'vza', vza, replace_nan=True)
                if verbosity > 1: print('Wrote vza')
                ac.output.nc_write(ofile, 'sza', sza, replace_nan=True)
                if verbosity > 1: print('Wrote sza')
                sza = None
                vza = None
        else:
            mus = np.asarray(gatts['mus'])  ## average cos sun zenith
            #mus.shape+=(1,1)

        ## write lat/lon
        if (output_geolocation):
            if os.path.exists(ofile):
                datasets = ac.shared.nc_datasets(ofile)
            else:
                datasets = []
            if ('lat' not in datasets) or ('lon' not in datasets):
                if verbosity > 1: print('Writing geolocation lon/lat')
                lon, lat = ac.shared.projection_geo(dct_prj)
                ac.output.nc_write(ofile, 'lon', lon, attributes=gatts, new=new, double=True)
                if verbosity > 1: print('Wrote lon')
                ac.output.nc_write(ofile, 'lat', lat, double=True)
                if verbosity > 1: print('Wrote lat')
                new=False

        ## write x/y
        if (output_xy):
            if os.path.exists(ofile):
                datasets = ac.shared.nc_datasets(ofile)
            else:
                datasets = []
            if ('x' not in datasets) or ('y' not in datasets):
                if verbosity > 1: print('Writing geolocation x/y')
                x, y = ac.shared.projection_geo(dct_prj, xy=True)
                ac.output.nc_write(ofile, 'x', x, new=new)
                if verbosity > 1: print('Wrote x')
                ac.output.nc_write(ofile, 'y', y)
                if verbosity > 1: print('Wrote y')
                new=False

        ## write TOA bands
        if verbosity > 1: print('Converting bands')
        for b in fmeta:
            if '.TIF' not in fmeta[b]['FILE']: continue
            if b in ['PIXEL', 'RADSAT']: continue
            if os.path.exists(fmeta[b]['FILE']):
                if b in waves_names:
                    pan = False
                    if b in pan_bands: ## pan band
                        if (not output_pan) & (not output_pan_ms): continue
                        pan = True
                        mus_pan = scipy.ndimage.zoom(mus, zoom=2, order=1) if len(np.atleast_1d(mus))>1 else mus * 1
                        data = ac.landsat.read_toa(fmeta[b], sub=sub_pan, mus=mus_pan, warp_to=warp_to_pan)
                        mus_pan = None
                    else: ## not a pan band
                        data = ac.landsat.read_toa(fmeta[b], sub=sub, mus=mus, warp_to=warp_to)
                    ds = 'rhot_{}'.format(waves_names[b])
                    ds_att = {'wavelength':waves_mu[b]*1000}
                    for k in fmeta[b]: ds_att[k] = fmeta[b][k]
                    if percentiles_compute:
                        ds_att['percentiles'] = percentiles
                        ds_att['percentiles_data'] = np.nanpercentile(data, percentiles)

                    if output_pan & pan:
                        ## write output
                        ofile_pan = ofile.replace('_L1R.nc', '_L1R_pan.nc')
                        ac.output.nc_write(ofile_pan, ds, data, attributes=gatts,replace_nan=True,
                                           new=new_pan, dataset_attributes = ds_att)
                        new_pan = False
                        if verbosity > 1: print('Converting bands: Wrote {} to separate L1R_pan'.format(ds))

                    ## prepare for low res output
                    if output_pan_ms & pan: data = scipy.ndimage.zoom(data, zoom=0.5, order=1)

                    ## write to ms file
                    ac.output.nc_write(ofile, ds, data, replace_nan=True, attributes=gatts, new=new, dataset_attributes = ds_att)
                    new = False
                    if verbosity > 1: print('Converting bands: Wrote {} ({})'.format(ds, data.shape))
                else:
                    if b in thermal_bands:
                        if output_thermal:
                            ds = 'bt{}'.format(b).lower()
                            ds_att = {'band':b}
                            for k in fmeta[b]: ds_att[k] = fmeta[b][k]
                            if percentiles_compute:
                                ds_att['percentiles'] = percentiles
                                ds_att['percentiles_data'] = np.nanpercentile(data, percentiles)
                            data = ac.landsat.read_toa(fmeta[b], sub=sub, warp_to=warp_to)
                            ac.output.nc_write(ofile, ds, data, replace_nan=True,
                                               attributes=gatts, new=new, dataset_attributes=ds_att)
                            new = False
                            if verbosity > 1: print('Converting bands: Wrote {}'.format(ds))
                    else:
                        continue

        if verbosity > 1:
            print('Conversion took {:.1f} seconds'.format(time.time()-t0))
            print('Created {}'.format(ofile))

        if limit is not None: sub = None
        if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles)
