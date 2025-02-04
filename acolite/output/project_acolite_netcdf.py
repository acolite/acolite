## def project_acolite_netcdf
## projects unprojected ACOLITE NetCDF data to a defined projection and extent
## written by Quinten Vanhellemont, RBINS
## 2022-01-04
## modifications: 2022-01-05 (QV) acolite function, changed handling provided x and y ranges
##                2022-01-10 (QV) renamed from reproject_acolite_netcdf
##                2022-07-05 (QV) determine projection limit from lat lon if none given
##                2022-07-06 (QV) simultaneous reprojection of multiple datasets (much faster!)
##                2022-10-17 (QV) added output_projection_polygon
##                2022-10-20 (QV) added nn option
##                2023-04-01 (QV) added viirs scanlines reprojection
##                2023-04-18 (QV) fixed reprojection for viirs overlapping  scanlines
##                2023-07-25 (QV) added counts output
##                2024-03-14 (QV) update settings handling
##                                added check for projection resolution
##                2024-04-16 (QV) update NetCDF writing, speed up read using gem
##                2025-01-30 (QV) use run instead of user settings
##                2025-02-04 (QV) improved settings handling

def project_acolite_netcdf(ncf, output = None, settings = None, target_file=None, output_counts = False):

    import os, time
    from pyproj import Proj

    import acolite as ac
    import numpy as np
    from pyresample.bilinear import NumpyBilinearResampler
    from pyresample import kd_tree, geometry

    ## read gem
    try:
        gem = ac.gem.gem(ncf)
        gem.store = False
    except:
        print('Error accessing {}, is this a NetCDF file?'.format(ncf))
        return

    gatts = gem.gatts
    if ('sensor' not in gatts):
        print('No sensor attribute in file {}'.format(ncf))
        return

    ## read datasets
    datasets = gem.datasets
    if ('lat' not in datasets) or ('lon' not in datasets):
        print('No lat/lon found in file {}'.format(ncf))
        return

    lon = None
    lat = None

    ## combine default and user defined settings
    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}
    ## get sensor specific defaults
    setd = ac.acolite.settings.parse(gem.gatts['sensor'])
    ## set sensor default if user has not specified the setting
    for k in setd:
        if k not in ac.settings['user']: setu[k] = setd[k]
    ## end set sensor specific defaults
    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    if (setu['output_projection_limit'] is None) & (setu['output_projection_polygon'] is not None):
        setu['output_projection_limit'] = ac.shared.polygon_limit(setu['output_projection_polygon'])

    if (setu['output_projection_limit'] is None) & ('limit_old' in setu):
        setu['output_projection_limit'] = [l for l in setu['limit_old']]

    if (setu['output_projection_limit'] is None) & (setu['limit'] is not None):
        setu['output_projection_limit'] = [l for l in setu['limit']]

    if (setu['output_projection_limit'] is None) & (setu['polygon'] is not None):
        setu['output_projection_limit'] = ac.shared.polygon_limit(setu['polygon'])

    if (setu['output_projection_limit'] is None):
        ## read lat/lon
        lat = gem.data('lat') #ac.shared.nc_data(ncf, 'lat')
        lon = gem.data('lat') #ac.shared.nc_data(ncf, 'lon')
        setu['output_projection_limit'] = [np.nanmin(lat), np.nanmin(lon), np.nanmax(lat), np.nanmax(lon)]
        print('Output limit from lat/lon', setu['output_projection_limit'])

    if setu['output_projection_resolution'] is not None:
        if len(setu['output_projection_resolution']) != 2:
            print('Provide a two element target_pixel_size.')
            return
        else:
            target_pixel_size = [float(v) for v in setu['output_projection_resolution']]
            target_pixel_size[1] *= -1
    else:
        if setu['default_projection_resolution'] is not None:
            target_pixel_size = [float(v) for v in setu['default_projection_resolution']]
            target_pixel_size[1] *= -1
            print('Using default grid size: {}x{}metres'.format(target_pixel_size[0], target_pixel_size[1]))
        else:
            print('No projection grid resolution set in output_projection_resolution or default_projection_resolution')
            return

    if setu['output_projection_limit'] is not None:
        if len(setu['output_projection_limit']) != 4:
            print('Provide a four element output_projection_limit.')
            return
        else:
            limit = [float(v) for v in setu['output_projection_limit']]

        if (setu['output_projection_epsg'] is None) & \
           (setu['output_projection_proj4'] is None):
              lon_ = (limit[1]+limit[3])/2
              lat_ = (limit[0]+limit[2])/2
              utm_zone, epsg = ac.shared.utm_epsg(lon_, lat_)
              print(utm_zone, epsg)
              setu['output_projection_epsg'] = epsg
    else:
        if not setu['output_projection_metres']:
            print('Provide a four element output_projection_limit.')
            return

    ## projection
    if setu['output_projection_epsg'] is not None:
        #projection = '+init=EPSG:{}'.format(setu['output_projection_epsg'])
        if 'EPSG' not in setu['output_projection_epsg']:
            projection = 'EPSG:{}'.format(setu['output_projection_epsg'])
        else:
            projection = '{}'.format(setu['output_projection_epsg'])
    elif setu['output_projection_proj4'] is not None:
        projection = setu['output_projection_proj4']
    else:
        print('No EPSG or proj4 string provided.')
        return

    ## user provided x and yrange
    if setu['output_projection_metres']:
        xrange_region = setu['output_projection_xrange']
        yrange_region = setu['output_projection_yrange']
        if (xrange_region is None) or  (yrange_region is None):
            print('Provide a output_projection_xrange and output_projection_yrange.')
            return
        if len(xrange_region) != 2:
            print('Provide a two element output_projection_xrange.')
            return
        if len(yrange_region) != 2:
            print('Provide a two element output_projection_yrange.')
            return

    ## create output file name
    bn = os.path.basename(ncf)
    bd = os.path.dirname(ncf)
    oname, nc = os.path.splitext(bn)
    if output is not None:
        bd = '{}'.format(output)
    elif setu['output'] is not None:
        bd = '{}'.format(setu['output'])
    ## add requested name or "reprojected"
    oname = '{}_{}'.format(oname, setu['output_projection_name'] if setu['output_projection_name'] is not None else "projected")
    ncfo = '{}/{}{}'.format(bd, oname, nc)

    print('Setting up target projection.')
    p = Proj(projection)

    ## find region extent
    if not setu['output_projection_metres']:
        ## project lat lon to metres
        xrange_raw, yrange_raw = p((limit[1],limit[1],limit[3],limit[3]),
                                   (limit[0],limit[2],limit[2],limit[0]))
        xrange_raw = (min(xrange_raw), max(xrange_raw))
        yrange_raw = (min(yrange_raw), max(yrange_raw))
        xrange_region = [xrange_raw[0] - (xrange_raw[0] % target_pixel_size[0]*2), xrange_raw[1]+target_pixel_size[0]*2-(xrange_raw[1] % target_pixel_size[0]*2)]
        yrange_region = [yrange_raw[1]+target_pixel_size[1]*2-(yrange_raw[1] % target_pixel_size[1]*2), yrange_raw[0] - (yrange_raw[0] % target_pixel_size[1]*2)]

    ## align grid to pixel size
    if setu['output_projection_resolution_align']:
        x_grid_off = xrange_region[0]%target_pixel_size[0], xrange_region[1]%target_pixel_size[0]
        y_grid_off = yrange_region[0]%target_pixel_size[1], yrange_region[1]%target_pixel_size[1]
        xrange = (xrange_region[0]-x_grid_off[0]), (xrange_region[1]+(target_pixel_size[0]-x_grid_off[1]))
        yrange = (yrange_region[0]-y_grid_off[0]), (yrange_region[1]+(target_pixel_size[0]-y_grid_off[1]))
    else:
        xrange = [xrange_region[0], xrange_region[1]]
        yrange = [yrange_region[0], yrange_region[1]]

    ## pixel sizes
    ny = int((yrange[1] - yrange[0])/target_pixel_size[1])
    nx = int((xrange[1] - xrange[0])/target_pixel_size[0])
    print(xrange, yrange)
    print(nx, ny)

    ## set up projection dict and nc_projection
    dct = {'xrange': xrange, 'yrange': yrange, 'p': p,
           'pixel_size': target_pixel_size, 'xdim': nx, 'ydim': ny}
    nc_projection = ac.shared.projection_netcdf(dct, add_half_pixel=True)

    ## set up target definition
    target_definition = geometry.AreaDefinition('area_id', 'description', 'proj_id',
                                          projection, nx, ny, [xrange[0],yrange[1],xrange[1],yrange[0]])

    ## read lat/lon
    if lat is None: lat = gem.data('lat') # ac.shared.nc_data(ncf, 'lat')
    if lon is None: lon = gem.data('lon') # ac.shared.nc_data(ncf, 'lon')

    ## make new output attributes
    gatts_out = {k:gatts[k] for k in gatts}
    for k in dct:
        if k == 'p': continue
        gatts_out[k] = dct[k]

    gatts_out['projection_key'] = [k for k in nc_projection if k not in ['x', 'y']][0]
    ## update oname in gatts
    gatts_out['oname'] = oname

    ## read in data to be reprojected
    data_in_stack = None
    datasets_out = []
    datasets_att = {}
    for ds in datasets:
        #data_in, att = ac.shared.nc_data(ncf, ds, attributes=True)
        data_in, att = gem.data(ds, attributes = True)
        if len(data_in.shape) != 2: continue
        if setu['verbosity'] > 2: print('Reading {} {}x{}'.format(ds, data_in.shape[0], data_in.shape[1]))
        if data_in_stack is None:
            data_in_stack = data_in
        else:
            data_in_stack = np.dstack((data_in_stack, data_in))
        datasets_out.append(ds)
        datasets_att[ds] = att
        data_in = None

    ## close input dataset
    gem.close()
    gem = None

    ## for NN/bilinear projection
    radius = setu['output_projection_radius']
    epsilon = setu['output_projection_epsilon']

    ## for bilinear projection
    neighbours = setu['output_projection_neighbours']

    if (setu['viirs_scanline_projection']) & ('VIIRS' in gatts['sensor']):
        slines = setu['viirs_scanline_width']
        if 'viirs_slines' in gatts: slines = gatts['viirs_slines']
        nscans = int(data_in_stack.shape[0]/slines)
        print('Assuming {} scans of {} lines'.format(nscans, slines))
        data_out_stack = np.zeros((ny,nx,len(datasets_out)))
        data_out_counts = np.zeros((ny,nx,len(datasets_out)), dtype=np.int8)
    else:
        nscans = 1

    ## reproject dataset
    t0 = time.time()
    print('Projecting datasets {}x{}x{} to {} {}x{}x{}'.format(data_in_stack.shape[0], data_in_stack.shape[1], len(datasets_out),\
                                                projection, nx, ny, len(datasets_out)))

    ## run through scans (1 if reprojecting all at once)
    for i in range(nscans):
        ## one reprojection
        if nscans == 1:
            ## set up source definition
            source_definition = geometry.SwathDefinition(lons=lon, lats=lat)
            ## set up resampler
            if setu['output_projection_resampling_method'] == 'bilinear':
                resampler = NumpyBilinearResampler(source_definition, target_definition, target_pixel_size[0]*radius,
                                                    neighbours=neighbours, epsilon=epsilon, reduce_data=False)
                data_out_stack = resampler.resample(data_in_stack, fill_value=np.nan)
            elif setu['output_projection_resampling_method'] == 'nearest':
                data_out_stack = kd_tree.resample_nearest(source_definition, data_in_stack, target_definition,
                                                            radius_of_influence=target_pixel_size[0]*radius,
                                                            epsilon=epsilon, fill_value=np.nan)
        ## reprojection per scan
        else:
            print('Reprojecting scan {}/{}'.format(i+1, nscans), end='\r')
            ts0 = time.time()
            ## set up scan source
            source_definition = geometry.SwathDefinition(lons=lon[i*slines:(i+1)*slines,:],
                                                         lats=lat[i*slines:(i+1)*slines,:])
            ## set up resampler
            if setu['output_projection_resampling_method'] == 'bilinear':
                resampler = NumpyBilinearResampler(source_definition, target_definition, target_pixel_size[0]*radius,
                                                    neighbours=neighbours, epsilon=epsilon, reduce_data=False)
                data_out_scan = resampler.resample(data_in_stack[i*slines:(i+1)*slines,:,:], fill_value=np.nan)
            elif setu['output_projection_resampling_method'] == 'nearest':
                data_out_scan = kd_tree.resample_nearest(source_definition, data_in_stack[i*slines:(i+1)*slines,:,:], target_definition,
                                                            radius_of_influence=target_pixel_size[0]*radius,
                                                            epsilon=epsilon, fill_value=np.nan)
            ## put scans in out stack
            scan_sub = np.where(np.isfinite(data_out_scan))
            data_out_stack[scan_sub] += data_out_scan[scan_sub]
            data_out_counts[scan_sub] += 1
            data_out_scan = None
            print('Reprojecting scan {}/{} took {:.1f} seconds'.format(i+1, nscans, time.time()-ts0), end='\r')
    if nscans >1: print()

    ## compute average overlapping area
    if nscans >1:
        data_out_stack/=data_out_counts
        data_out_stack[data_out_counts==0] = np.nan
        if not output_counts: data_out_counts = None

    data_in_stack = None
    if setu['output_projection_fillnans']:
        data_out_stack[data_out_stack == 0] = np.nan
        if len(data_out_stack.shape) == 3:
            for di in range(data_out_stack.shape[2]):
                data_out_stack[:,:,di] = ac.shared.fillnan(data_out_stack[:,:,di], max_distance=setu['output_projection_filldistance'])
        else:
            data_out_stack[:,:] = ac.shared.fillnan(data_out_stack[:,:], max_distance=setu['output_projection_filldistance'])
    t1 = time.time()
    print('Reprojection of {} datasets took {:.1f} seconds'.format(len(datasets_out),t1-t0))

    ## setup output dataset
    gemo = ac.gem.gem(ncfo, new = True)
    gemo.gatts = {k: gatts_out[k] for k in gatts_out} ## set output attributes
    gemo.nc_projection = nc_projection ## set new projection details
    gemo.verbosity = 6

    ## write results
    for di, ds in enumerate(datasets_out):
        if ds in ['l2_flags']: continue

        if setu['verbosity'] > 2: print('Writing {} {}x{}'.format(ds, data_out_stack[:,:,di].shape[0], data_out_stack[:,:,di].shape[1]))
        gemo.write(ds, data_out_stack[:,:,di], ds_att = datasets_att[ds])
        if output_counts:
            if setu['verbosity'] > 2: print('Writing {} {}x{}'.format(ds+'_n', data_out_counts[:,:,di].shape[0], data_out_counts[:,:,di].shape[1]))
            gemo.write(ds+'_n', data_out_counts[:,:,di])
    data_out_stack = None
    data_out_counts = None

    ## compute target lon/lat
    tlon, tlat = ac.shared.projection_geo(dct)
    gemo.write('lon', tlon)
    tlon = None
    gemo.write('lat', tlat)
    tlat = None

    ## close output dataset
    gemo.close()
    gemo = None

    print('Wrote {}'.format(ncfo))
    return(ncfo)
