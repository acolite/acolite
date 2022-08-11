## def agh
## ACOLITE/GEE Hybrid processing
## extracts either L1R data from GEE or performs modified fixed DSF and gets L2R
## written by Quinten Vanhellemont, RBINS
## 2022-04-11
## modifications: 2022-04-12 (QV) added NetCDF output, added resolved geometry output
##                2022-04-13 (QV) added to acolite, new default naming of gems, removed limit and limit_name keywords
##                                added tiled download of large files
##                2022-04-14 (QV) added glint correction
##                2022-04-15 (QV) added ancillary data
##                2022-05-22 (QV) added download of BT data from Landsat, added metadata copy
##                2022-07-18 (QV) added check for existing files & override setting

def agh(image, imColl, rsrd = {}, lutd = {}, luti = {}, settings = {}):
    import os, dateutil.parser, requests
    import acolite as ac
    from acolite import gee ## currently not imported in main acolite

    import numpy as np
    from osgeo import gdal
    gdal.UseExceptions()

    import ee
    #ee.Authenticate() ## assume ee use is authenticated in current environment
    ee.Initialize()

    uoz = settings['uoz_default']
    uwv = settings['uwv_default']
    pressure = settings['pressure']
    wind = settings['wind']

    rhop_par = settings['rhop_par']
    sel_par = settings['sel_par']

    limit = None
    region_name = None
    output = os.getcwd()
    if 'limit' in settings: limit = settings['limit']
    if 'region_name' in settings: region_name = settings['region_name']
    if 'output' in settings: output = settings['output']
    if not os.path.exists(output): os.makedirs(output)

    if limit is not None:
        rname = 'crop' if region_name is None else region_name
    else:
        rname = ''

    fkey, pid = image[0], image[1]

    tar_band = 'B2' ## target band for projection and output resolution 10m for S2, 30m for Landsat
    if ('LANDSAT' in fkey) & (pid[0:3] == 'LT0'): tar_band = 'B10' ## TIRS only data

    ## output file names for combined file
    ## used to determine NetCDF name, also for Google Drive outputs
    ext = ''
    if len(rname) >= 0: ext = '_{}'.format(rname)
    rhot_file = pid+'_rhot'+ext
    rhot_file_local = '{}/{}.zip'.format(output,rhot_file)
    rhos_file = pid+'_rhos'+ext
    rhos_file_local = '{}/{}.zip'.format(output,rhos_file)
    geom_file = pid+'_geom'+ext
    geom_file_local = '{}/{}.zip'.format(output,geom_file)

    file_type = 'L1R_GEE'
    if settings['run_hybrid_dsf']: file_type = 'L2R_GEE_HYBRID'

    ## select product
    i = imColl.filter(ee.Filter.eq(fkey, pid)).first()

    ## get projection info
    ## get projection
    proj = i.select(tar_band).projection().getInfo()
    if 'crs' in proj:
        proj_crs = proj['crs']
    elif 'wkt' in proj:
        proj_crs = proj['wkt']
    p = ee.Projection(proj_crs, proj['transform'])
    scale = p.nominalScale().getInfo()
    if settings['output_scale'] is not None: scale = settings['output_scale']

    if limit is not None:
        if settings['strict_subset']:
            ## determin strict lat/lon rectangle box
            region = ee.Geometry.BBox(limit[1], limit[0], limit[3], limit[2])
        else:
            ## determine image bounding box
            imx = []
            imy = []
            for ii in [[1,0], [1,2], [3,0], [3,2]]:
                ## make point geometry
                pt = ee.Geometry.Point([limit[ii[0]], limit[ii[1]]])

                ## get pixel coordinates in x/y
                tmp = ee.Image.clip(i.pixelCoordinates(i.select(tar_band).projection()),pt)
                ret = ee.Image.reduceRegion(tmp, ee.Reducer.toList()).getInfo()
                imx.append(ret['x'][0])
                imy.append(ret['y'][0])

            ## use pixel coordinates from image to make new subset
            eesub = ee.List([min(imx), min(imy), max(imx), max(imy)])
            region = ee.Geometry.Rectangle(eesub, p, True, False)

    ## subset here if local aot is to be computed
    if (limit is not None) & (settings['subset_aot']): i = i.clip(region)

    im = i.getInfo()
    bands = [p['id'] for p in im['bands'] if p['id'][0]]

    if 'PRODUCT_ID' in im['properties']: ## Sentinel-2 image
        fkey = 'PRODUCT_ID'
        pid = im['properties'][fkey]
        dtime = [k for k in im['properties']['GRANULE_ID'].split('_') if len(k) == 15][0]
        tile_name = '{}'.format(im['properties']['MGRS_TILE'])
        satellite = pid[0:3]
        sensor = '{}_MSI'.format(satellite)
        scale_factor = 0.0001
        add_factor = 0
        if (im['properties']['PROCESSING_BASELINE'][1]>='4') & ~('S2_HARMONIZED' in im['id']):
            add_factor = -1000
    elif 'LANDSAT_PRODUCT_ID' in im['properties']: ## Landsat image
        fkey = 'LANDSAT_PRODUCT_ID'
        pid = im['properties'][fkey]
        dtime = im['properties']['DATE_ACQUIRED']+'T'+im['properties']['SCENE_CENTER_TIME']
        row = '{}'.format(im['properties']['WRS_ROW']).zfill(3)
        path = '{}'.format(im['properties']['WRS_PATH']).zfill(3)
        tile_name = '{}{}'.format(path, row)
        satellite = pid[0]+pid[3]
        if satellite == 'L5':
            sensor = '{}_TM'.format(satellite)
        elif satellite == 'L7':
            sensor = '{}_ETM'.format(satellite)
        elif satellite in ['L8', 'L9']:
            sensor = '{}_OLI'.format(satellite)
        scale_factor = 1
        add_factor = 0

    ## parse datetime
    dt = dateutil.parser.parse(dtime)

    if settings['use_scene_name']:
        ofile = rhot_file_local.replace('_rhot', '_{}'.format(file_type)).replace('.zip', '.nc')
    else:
        obase  = '{}_{}_{}_{}{}'.format(sensor,  dt.strftime('%Y_%m_%d_%H_%M_%S'), tile_name, file_type, ext)
        ofile = '{}/{}.nc'.format(os.path.dirname(rhot_file_local), obase)

    if os.path.exists(ofile) & (settings['override'] is False):
        print('AGH file {} exists, set override=True to replace'.format(ofile))
        return()

    ## get ancillary
    if settings['ancillary_data']:
        print('Getting ancillary for scene centre location.')
        ## get central lon and lat
        mp = i.pixelLonLat().reproject(p, None, scale * 1.0)
        ll = mp.reduceRegion(geometry=(i.geometry() if limit is None else region), \
                             reducer= ee.Reducer.mean(), bestEffort=True).getInfo()
        print('Scene centre: {:.5f}E {:.5f}N'.format(ll['longitude'], ll['latitude']))
        anc = ac.ac.ancillary.get(dt, ll['longitude'], ll['latitude'])

        ## overwrite the defaults
        if ('ozone' in anc): uoz = anc['ozone']['interp']/1000. ## convert from MET data
        if ('p_water' in anc): uwv = anc['p_water']['interp']/10. ## convert from MET data
        if ('z_wind' in anc) & ('m_wind' in anc) & (wind is None):
            wind = ((anc['z_wind']['interp']**2) + (anc['m_wind']['interp']**2))**0.5
        if ('press' in anc): pressure = anc['press']['interp']
        print(uoz, uwv, wind, pressure)

    if sensor not in rsrd:
        print('Loading RSRs: {}'.format(sensor))
        rsrd[sensor] = ac.shared.rsr_dict(sensor=sensor)[sensor]

    if (sensor not in lutd) & (settings['run_hybrid_dsf']):
        print('Loading atmosphere LUTs: {}'.format(sensor))
        lutd[sensor] = ac.aerlut.import_luts(sensor=sensor)

    if (sensor not in luti) & (settings['run_hybrid_dsf']) & (settings['glint_correction']):
        print('Loading interface LUTs: {}'.format(sensor))
        luti[sensor] = ac.aerlut.import_rsky_luts(models=[1,2], lutbase='ACOLITE-RSKY-202102-82W', sensor=sensor)

    ## get average geometry
    ## will be updated later if SAA/SZA/VAA/VZA data are available
    geometry = {}
    if sensor in ['S2A_MSI', 'S2B_MSI']:
        geometry['saa'] = im['properties']['MEAN_SOLAR_AZIMUTH_ANGLE']
        geometry['sza'] = im['properties']['MEAN_SOLAR_ZENITH_ANGLE']
        bgeometry = {}
        for b in bands:
            if 'QA' in b: continue
            if 'SCL' in b: continue
            if 'MEAN_INCIDENCE_AZIMUTH_ANGLE_{}'.format(b) in im['properties']:
                vaa = im['properties']['MEAN_INCIDENCE_AZIMUTH_ANGLE_{}'.format(b)]
            else: continue
            if 'MEAN_INCIDENCE_ZENITH_ANGLE_{}'.format(b) in im['properties']:
                vza = im['properties']['MEAN_INCIDENCE_ZENITH_ANGLE_{}'.format(b)]
            else: continue
            bgeometry[b] = {'vza':vza, 'vaa':vaa}
        geometry['vza'] = np.nanmean([bgeometry[b]['vza'] for b in bgeometry])
        geometry['vaa'] = np.nanmean([bgeometry[b]['vaa'] for b in bgeometry])
        aot_skip_bands = ['9', '10']
        glint_bands = ['11', '12']
    else:
        geometry['saa'] = im['properties']['SUN_AZIMUTH']
        geometry['sza'] = 90 - im['properties']['SUN_ELEVATION']
        geometry['vza'] = 0
        geometry['vaa'] = 0
        aot_skip_bands = ['9']
        if sensor in ['L8_OLI', 'L9_OLI']:
            glint_bands = ['6', '7']
        else:
            glint_bands = ['5', '7']

    ## make rhot dataset
    i_rhot = None
    for b in rsrd[sensor]['rsr_bands']:
        bname = 'B{}'.format(b)
        print('rhot {}'.format(bname))
        ## gas and path corrected rhot
        rhot = i.expression('(DN + {}) * {}'.format(add_factor, scale_factor), {'DN': i.select(bname)})
        if i_rhot is None:
            i_rhot = ee.Image(rhot)
        else:
            i_rhot = i_rhot.addBands(rhot)

    ## add thermal bands
    ## not rhot but at sensor BT
    thermal_bands = []
    for b in im['bands']:
        bname = b['id']
        if (bname[0] == 'B') & (bname[1:] not in rsrd[sensor]['rsr_bands']):
            print(bname)
            bt = i.expression('(DN + {}) * {}'.format(add_factor, scale_factor), {'DN': i.select(bname)})
            i_rhot = i_rhot.addBands(bt)
            thermal_bands.append(bname)

    ## TOA band info
    rhot_info = i_rhot.getInfo()
    obands_rhot = [ib['id'] for ib in rhot_info['bands']]

    ## make geometry per pixel dataset if available
    i_geom = None
    for b in ['SAA', 'SZA', 'VAA', 'VZA']:
        if b not in bands: continue
        print(b)
        ang = i.expression('DN/100', {'DN': i.select(b)})
        if i_geom is None:
            i_geom = ee.Image(ang)
        else:
            i_geom = i_geom.addBands(ang)
    if i_geom is not None:
        geom_info = i_geom.getInfo()
        obands_geom = [ib['id'] for ib in geom_info['bands']]

    print('Getting geometry and gas transmittance')
    ## get geometry percentiles
    geom_percentile = 50
    prc = i.reduceRegion(reducer= ee.Reducer.percentile([geom_percentile]), bestEffort=True).getInfo()
    for b in bands:
        ### get geometry if available
        if b in ['SAA', 'SZA', 'VAA', 'VZA']:
            k = '{}_p{}'.format(b,geom_percentile)
            if k in prc: geometry[b.lower()] = prc[k]/100

    ## update geometry
    geometry['raa'] = np.abs(geometry['saa'] - geometry['vaa'])
    while geometry['raa'] > 180:
        geometry['raa'] = 360 - geometry['raa']

    ## get gas transmittance
    ttg = ac.ac.gas_transmittance(geometry['sza'], geometry['vza'], pressure = pressure,
                            uoz = uoz, uwv = uwv, rsr=rsrd[sensor]['rsr'])

    ## run fixed DSF on this data
    if settings['run_hybrid_dsf']:
        ## get rhot percentiles
        percentiles = settings['percentiles']
        print('Extracting rhot percentiles')
        prc = i_rhot.reduceRegion(reducer= ee.Reducer.percentile(percentiles), bestEffort=True).getInfo()
        prc_data = {p: {b: prc['{}_p{}'.format(b,p)] for b in obands_rhot} for p in percentiles}

        ## fit aerosol models
        print('Fitting aerosol models')
        results = {}
        luts = list(lutd[sensor].keys())
        for lut in luts:
            taua_arr = None
            rhot_arr = None
            taua_bands = []
            ## run through bands
            for b in rsrd[sensor]['rsr_bands']:
                if b in aot_skip_bands: continue
                print(sensor, b, lut)
                ret = lutd[sensor][lut]['rgi'][b]((pressure, lutd[sensor][lut]['ipd']['romix'],
                                           geometry['raa'], geometry['vza'], geometry['sza'], lutd[sensor][lut]['meta']['tau']))
                rhot = np.asarray([prc_data[p]['B{}'.format(b)] for p in percentiles])
                rhot /= ttg['tt_gas'][b]
                taua = np.interp(rhot, ret, lutd[sensor][lut]['meta']['tau'])
                if taua_arr is None:
                    rhot_arr = 1.0 * rhot
                    taua_arr = 1.0 * taua
                else:
                    rhot_arr = np.vstack((rhot_arr, rhot))
                    taua_arr = np.vstack((taua_arr, taua))

                taua_bands.append(b)

            ## find aot value
            bidx = np.argsort(taua_arr[:,settings['pidx']])
            taua = np.nanmean(taua_arr[bidx[0:settings['nbands']],settings['pidx']])
            taua_std = np.nanstd(taua_arr[bidx[0:settings['nbands']],settings['pidx']])
            taua_cv = taua_std/taua
            taua, taua_std, taua_cv*100

            ## store results
            results[lut] = {'taua_bands': taua_bands, 'taua_arr': taua_arr, 'rhot_arr': rhot_arr,
                            'taua': taua, 'taua_std': taua_std, 'taua_cv': taua_cv,'bidx': bidx}

        ## select LUT and aot
        sel_lut = None
        sel_aot = None
        sel_val = np.inf
        for lut in results:
            if results[lut][sel_par] < sel_val:
                sel_val = results[lut][sel_par] * 1.0
                sel_aot = results[lut]['taua'] * 1.0
                sel_lut = '{}'.format(lut)

        ## get atmosphere parameters
        print('Getting final atmosphere parameters')
        am = {}
        for par in lutd[sensor][sel_lut]['ipd']:
            am[par] = {b: lutd[sensor][sel_lut]['rgi'][b]((pressure, lutd[sensor][sel_lut]['ipd'][par],
                                                           geometry['raa'], geometry['vza'], geometry['sza'], sel_aot))
                       for b in rsrd[sensor]['rsr_bands']}

        print(sel_lut, sel_aot, sel_val)

        ## do atmospheric correction
        print('Performing atmospheric correction')
        i_rhos = None
        for b in rsrd[sensor]['rsr_bands']:
            #if b in aot_skip_bands: continue
            bname = 'B{}'.format(b)

            ## gas and path corrected rhot
            rhot_prime = i_rhot.expression('(rhot / {}) - {}'.format(ttg['tt_gas'][b], am[rhop_par][b]), \
                                           {'rhot': i_rhot.select(bname)})

            ## rhos
            rhos = rhot_prime.expression('(rhotp) / ({} - {} * rhotp)'.format(am['dutott'][b], am['astot'][b]), \
                                         {'rhotp': rhot_prime.select(bname)})

            if i_rhos is None:
                i_rhos = ee.Image(rhos)
            else:
                i_rhos = i_rhos.addBands(rhos)

        ## do atmospheric correction + glint
        if settings['glint_correction']:
            print('Performing glint correction using bands {} {}'.format(glint_bands[0], glint_bands[1]))
            model = int(sel_lut[-1])
            g1 = i_rhos.select('B{}'.format(glint_bands[0]))
            g2 = i_rhos.select('B{}'.format(glint_bands[1]))
            glint = (g1.add(g2)).divide(2).rename('glint')
            glintMask = glint.gt(settings['glint_min']).multiply(glint.lt(settings['glint_max'])).selfMask()
            #glint = glint.multiply(glintMask)
            glint = glint.mask(glintMask)
            glint = glint.unmask(0)
            glint_dict = {b:luti[sensor][model]['rgi'][b]((geometry['raa'], geometry['vza'], geometry['sza'], settings['glint_wind'], sel_aot)) for b in rsrd[sensor]['rsr_bands']}
            glint_spec = np.asarray([glint_dict[b] for b in glint_dict])
            glint_ave = {b: glint_dict[b]/((glint_dict[glint_bands[0]]+glint_dict[glint_bands[1]])/2) for b in glint_dict}
            i_rhos = None
            for b in rsrd[sensor]['rsr_bands']:
                #if b in aot_skip_bands: continue
                bname = 'B{}'.format(b)

                ## gas and path corrected rhot
                rhot_prime = i_rhot.expression('(rhot / {}) - {}'.format(ttg['tt_gas'][b], am[rhop_par][b]), \
                                               {'rhot': i_rhot.select(bname)})

                ## rhos not corrected for glint
                rhos = rhot_prime.expression('(rhotp) / ({} - {} * rhotp)'.format(am['dutott'][b], am['astot'][b]), \
                                             {'rhotp': rhot_prime.select(bname)})

                ## rhos corrected for glint
                rhos = rhos.expression('rhos - (rhog * {})'.format(glint_ave[b]), \
                                             {'rhos': rhos.select(bname), 'rhog': glint.select('glint')})

                if i_rhos is None:
                    i_rhos = ee.Image(rhos)
                else:
                    i_rhos = i_rhos.addBands(rhos)
        rhos_info = i_rhos.getInfo()
        obands_rhos = [ib['id'] for ib in rhos_info['bands']]

    ## find image dimensions
    for b in rhot_info['bands']:
        if b['id'] == tar_band:
            odim = b['dimensions']
            origin = b['origin']
            print('Region dimensions {}'.format(odim))
            print('Region origin {}'.format(origin))

    ## check size and tile if necessary
    tiled_transfer = False
    tile_size = [606,606]
    tiles_grid = 1,1
    if ((odim[0]*odim[1]*13)*9 > 50331648) or (odim[0]>tile_size[0]) or (odim[1]>tile_size[1]): #((odim[0]*odim[1]) > (tile_size[0]*tile_size[1])):
        print('Warning large dataset: {}x{}'.format(odim[0],odim[1]))
        tiled_transfer = True
        tiles_grid = int(np.ceil(odim[0]/tile_size[0])), int(np.ceil(odim[1]/tile_size[1]))
    num_tiles = tiles_grid[0]*tiles_grid[1]

    ## identify tiles
    if tiled_transfer:
        tiles = []
        tn = 0
        for ti in range(tiles_grid[0]):
            for tj in range(tiles_grid[1]):
                ti0 = ti*tile_size[0]
                ti1 = min((ti+1)*tile_size[0], odim[0])
                tj0 = tj*tile_size[1]
                tj1 = min((tj+1)*tile_size[1], odim[1])
                tiles.append([ti0, ti1, tj0, tj1, 't{}'.format(str(tn).zfill(3))])
                print('tile', tiles[-1])
                tn+=1
    else:
        tiles = [[0,odim[0], 0,odim[1],'']]
    print('Running download with {} tiles'.format(len(tiles)))

#    ## output file names
#    ext = ''
#    if len(rname) >= 0: ext = '_{}'.format(rname)
#    rhot_file = pid+'_rhot'+ext
#    rhot_file_local = '{}/{}.zip'.format(output,rhot_file)
#    rhos_file = pid+'_rhos'+ext
#    rhos_file_local = '{}/{}.zip'.format(output,rhos_file)
#    geom_file = pid+'_geom'+ext
#    geom_file_local = '{}/{}.zip'.format(output,geom_file)

    ## set output config for rhot
    output_config = {'description': rhot_file, 'scale': scale,  'folder':'ACOLITE'}

    ## test if setting crs works
    output_config['crs'] = proj_crs
    output_config['crs_transform'] = proj['transform']

    ## set output dir
    if settings['drive_output'] is not None: output_config['folder'] = settings['drive_output']
    if limit is not None:
        output_config['region'] = region
    else:
        output_config['region'] = i.geometry()

    ## output data to drive
    if settings['store_output_google_drive']:
        if settings['store_rhot']:
            output_config['description'] = rhot_file
            task_rhot = ee.batch.Export.image.toDrive(i_rhot.select(obands_rhot), **output_config)
            print('Exporting to Drive {}'.format(output_config['description']))
            task_rhot.start()
            status = gee.check_task(task_rhot, task_sleep = settings['task_check_sleep'])

        if settings['store_geom'] & (i_geom is not None):
            output_config['description'] = geom_file
            task_geom = ee.batch.Export.image.toDrive(i_geom.select(obands_geom), **output_config)
            print('Exporting to Drive {}'.format(output_config['description']))
            task_geom.start()
            status = gee.check_task(task_geom, task_sleep = settings['task_check_sleep'])

        if settings['store_rhos'] & settings['run_hybrid_dsf']:
            output_config['description'] = rhos_file
            task_rhos = ee.batch.Export.image.toDrive(i_rhos.select(obands_rhos), **output_config)
            print('Exporting to Drive {}'.format(output_config['description']))
            task_rhos.start()
            status = gee.check_task(task_rhos, task_sleep = settings['task_check_sleep'])
    ## end output to drive

    ## store data locally, with tiling if needed
    if settings['store_output_locally']:
        ## empty list to track local files
        rhot_files = []
        geom_files = []
        rhos_files = []
        for ti, tile in enumerate(tiles):
            ## output file names
            ext = ''
            if len(rname) >= 0: ext = '_{}'.format(rname)
            if tiled_transfer:
                tile_id = tile[4]
                ext+='_'+tile_id
                print('Getting tile {} {}/{}'.format(tile_id, ti+1, num_tiles))
                ## use pixel coordinates to make tile subset
                mins = ee.List([origin[0]+tile[0], origin[1]+tile[2]])
                maxs = ee.List([origin[0]+tile[1], origin[1]+tile[3]])
                ## set output region to this tile
                rect = ee.Geometry.Rectangle(mins.cat(maxs), p, True, False)#.transform("EPSG:4326")
                output_config['region'] = rect

            ## set output names for this tile
            rhot_file_tile = pid+'_rhot'+ext
            rhot_file_tile_local = '{}/{}.zip'.format(output,rhot_file_tile)
            rhos_file_tile = pid+'_rhos'+ext
            rhos_file_tile_local = '{}/{}.zip'.format(output,rhos_file_tile)
            geom_file_tile = pid+'_geom'+ext
            geom_file_tile_local = '{}/{}.zip'.format(output,geom_file_tile)

            ## write rhot (tile)
            if settings['store_rhot']:
                output_config['description'] = rhot_file_tile
                url = i_rhot.getDownloadUrl({
                        'name': output_config['description'],
                        'bands': obands_rhot,
                        'region': output_config['region'],
                        'scale': output_config['scale'],
                        'crs': output_config['crs'],
                        'crs_transform': output_config['crs_transform'],
                        'filePerBand': False})
                print('Downloading {}'.format(rhot_file_tile))
                response = requests.get(url)
                with open(rhot_file_tile_local, 'wb') as f:
                    f.write(response.content)
                rhot_files.append(rhot_file_tile_local)

            ## write geometry (tile)
            if settings['store_geom'] & (i_geom is not None):
                output_config['description'] = geom_file_tile
                url = i_geom.getDownloadUrl({
                    'name': output_config['description'],
                    'bands': obands_geom,
                    'region': output_config['region'],
                    'scale': output_config['scale'],
                    'crs': output_config['crs'],
                    'crs_transform': output_config['crs_transform'],
                    'filePerBand': False})
                print('Downloading {}'.format(geom_file_tile))
                response = requests.get(url)
                with open(geom_file_tile_local, 'wb') as f:
                    f.write(response.content)
                geom_files.append(geom_file_tile_local)

            ## write rhos (tile)
            if settings['store_rhos'] & settings['run_hybrid_dsf']:
                output_config['description'] = rhos_file_tile
                url = i_rhos.getDownloadUrl({
                    'name': output_config['description'],
                    'bands': obands_rhos,
                    'region': output_config['region'],
                    'scale': output_config['scale'],
                    'crs': output_config['crs'],
                    'crs_transform': output_config['crs_transform'],
                    'filePerBand': False})
                print('Downloading {}'.format(rhos_file_tile))
                response = requests.get(url)
                with open(rhos_file_tile_local, 'wb') as f:
                    f.write(response.content)
                rhos_files.append(rhos_file_tile_local)
    ## end store local files

    ## output to ACOLITE style NetCDF
    if settings['convert_output']:
        verbosity = 5
        min_tgas_rho = 0.75
        new = True
        setu = {'netcdf_compression': True, 'netcdf_compression_level': 4}

        ## read data
        rhot_data = None
        rhos_data = None
        geom_data = None
        for ti, tile in enumerate(tiles):
            print(tile[4])
            rhotf = rhot_files[ti]
            rhot_image_file = '/vsizip/{}/{}.tif'.format(rhotf, os.path.basename(os.path.splitext(rhotf)[0]))
            ## read rhot
            if os.path.exists(rhotf):
                if ti == 0:
                    dct = ac.shared.projection_read(rhot_image_file)
                else:
                    ## update dct
                    dct_ = ac.shared.projection_read(rhot_image_file)
                    dct['yrange'] = max(dct['yrange'][0], dct_['yrange'][0]), min(dct['yrange'][1], dct_['yrange'][1])
                    dct['xrange'] = min(dct['xrange'][0], dct_['xrange'][0]), max(dct['xrange'][1], dct_['xrange'][1])
                    dct['ydim'] = int((dct['yrange'][1]-dct['yrange'][0])/dct['pixel_size'][1])
                    dct['xdim'] = int((dct['xrange'][1]-dct['xrange'][0])/dct['pixel_size'][0])
                    dct['dimensions'] = (dct['xdim'], dct['ydim'])

                ds = gdal.Open(rhot_image_file)
                if num_tiles == 1:
                    rhot_data = ds.ReadAsArray()
                else:
                    if rhot_data is None: rhot_data = np.zeros((ds.RasterCount, odim[1], odim[0]))+np.nan
                    rhot_data[:, tile[2]:tile[3], tile[0]:tile[1]] = ds.ReadAsArray()
                ds = None

            ## read geom data
            if len(geom_files) == num_tiles:
                geomf = geom_files[ti]
                geom_image_file = '/vsizip/{}/{}.tif'.format(geomf, os.path.basename(os.path.splitext(geomf)[0]))
                ## read geom
                if os.path.exists(geomf):
                    ds = gdal.Open(geom_image_file)
                    if num_tiles == 1:
                        geom_data = ds.ReadAsArray()
                    else:
                        if geom_data is None: geom_data = np.zeros((ds.RasterCount, odim[1], odim[0]))+np.nan
                        geom_data[:, tile[2]:tile[3],tile[0]:tile[1]] = ds.ReadAsArray()
                    ds = None

            ## read rhos data
            if len(rhos_files) == num_tiles:
                rhosf = rhos_files[ti]
                rhos_image_file = '/vsizip/{}/{}.tif'.format(rhosf, os.path.basename(os.path.splitext(rhosf)[0]))
                ## read rhos
                if os.path.exists(rhosf):
                    ds = gdal.Open(rhos_image_file)
                    if num_tiles == 1:
                        rhos_data = ds.ReadAsArray()
                    else:
                        if rhos_data is None: rhos_data = np.zeros((ds.RasterCount, odim[1], odim[0]))+np.nan
                        rhos_data[:, tile[2]:tile[3],tile[0]:tile[1]] = ds.ReadAsArray()
                    ds = None


        ## image file in zip
        #rhot_image_file = '/vsizip/{}/{}.tif'.format(rhot_file_local, rhot_file)
        #rhos_image_file = '/vsizip/{}/{}.tif'.format(rhos_file_local, rhos_file)
        #geom_image_file = '/vsizip/{}/{}.tif'.format(geom_file_local, geom_file)

        ## gatts for output NetCDF
        gatts = {}
        gatts['acolite_type'] = 'l1r_gee'
        gatts['pid'] = pid
        gatts['sensor'] = sensor
        gatts['isodate'] = dt.isoformat()
        ## add image metadata from properties
        property_keys = list(im['properties'].keys())
        property_keys.sort()
        for k in property_keys: gatts[k] = im['properties'][k]
        ## add geometry average
        for k in geometry: gatts[k] = geometry[k]
        #file_type = 'L1R_GEE'

        if settings['run_hybrid_dsf']:
            gatts['acolite_type'] = 'l2r_gee_hybrid'
            gatts['pressure'] = pressure
            gatts['uwv'] = uwv
            gatts['uoz'] = uoz
            gatts['ac_aot_550'] = sel_aot
            gatts['ac_model'] = sel_lut
            gatts['ac_var'] = sel_val
            #file_type = 'L2R_GEE_HYBRID'

        ## output file names
        ext = ''
        if len(rname) >= 0: ext = '_{}'.format(rname)
        #rhot_file = pid+'_rhot'+ext
        #rhot_file_local = '{}/{}.zip'.format(output,rhot_file)
        if settings['use_scene_name']:
            ofile = rhot_file_local.replace('_rhot', '_{}'.format(file_type)).replace('.zip', '.nc')
        else:
            obase  = '{}_{}_{}_{}{}'.format(gatts['sensor'],  dt.strftime('%Y_%m_%d_%H_%M_%S'), tile_name, file_type, ext)
            ofile = '{}/{}.nc'.format(os.path.dirname(rhot_file_local), obase)
        print(ofile)

        ## convert dct info into nc_projection and lat, lon
        #dct = ac.shared.projection_read(rhot_image_file)
        nc_projection = ac.shared.projection_netcdf(dct, add_half_pixel=True)
        lon, lat = ac.shared.projection_geo(dct, add_half_pixel=True)

        ## write lon
        ac.output.nc_write(ofile, 'lon', lon, attributes=gatts, new=new, double=True, nc_projection=nc_projection,
                           netcdf_compression=setu['netcdf_compression'],
                           netcdf_compression_level=setu['netcdf_compression_level'])
        if verbosity > 1: print('Wrote lon')
        lon = None

        ## write lat
        ac.output.nc_write(ofile, 'lat', lat, double=True,
                           netcdf_compression=setu['netcdf_compression'],
                           netcdf_compression_level=setu['netcdf_compression_level'])
        if verbosity > 1: print('Wrote lat')
        lat = None
        new = False

        if geom_data is not None:
            for bi, b in enumerate(['SAA', 'SZA', 'VAA', 'VZA']):
                ac.output.nc_write(ofile, b.lower(), geom_data[bi,:,:],
                                         netcdf_compression=setu['netcdf_compression'],
                                         netcdf_compression_level=setu['netcdf_compression_level'])
                if verbosity > 1: print('Wrote {}'.format(b.lower()))
            ## compute raa
            raa = np.abs(geom_data[0,:,:]-geom_data[2,:,:])
            raa[raa>180] = np.abs(360 - raa[raa>180])
            ac.output.nc_write(ofile, 'raa', raa,
                                netcdf_compression=setu['netcdf_compression'],
                                netcdf_compression_level=setu['netcdf_compression_level'])
            if verbosity > 1: print('Wrote raa')
            raa = None
            geom_data = None

        for bi, b in enumerate(rsrd[sensor]['rsr_bands']):
            wave_name = rsrd[sensor]['wave_name'][b]
            wave_nm = rsrd[sensor]['wave_nm'][b]
            att = {'band': b, 'wave_name': wave_name, 'wave_nm': wave_nm}
            for k in ttg: att[k] = ttg[k][b]

            rhot_ds = 'rhot_{}'.format(wave_name)
            rhos_ds = 'rhos_{}'.format(wave_name)

            ## write rhot data
            if rhot_data is not None:
                ac.output.nc_write(ofile, rhot_ds, rhot_data[bi, :, :], dataset_attributes=att,
                                   netcdf_compression=setu['netcdf_compression'], netcdf_compression_level=setu['netcdf_compression_level'])
                print('Wrote {}'.format(rhot_ds))

            ## write rhos data
            if (att['tt_gas'] > min_tgas_rho) & (rhos_data is not None):
                ac.output.nc_write(ofile, rhos_ds, rhos_data[bi, :, :], dataset_attributes=att,
                                   netcdf_compression=setu['netcdf_compression'], netcdf_compression_level=setu['netcdf_compression_level'])
                print('Wrote {}'.format(rhos_ds))

        ## write thermal bands
        for bi, b in enumerate(obands_rhot):
            if b not in thermal_bands: continue
            dso = b.replace('B', 'bt')
            ac.output.nc_write(ofile, dso, rhot_data[bi, :, :], dataset_attributes={'band': b},
                               netcdf_compression=setu['netcdf_compression'], netcdf_compression_level=setu['netcdf_compression_level'])
            print('Wrote {}'.format(dso))


        print('Wrote {}'.format(ofile))
        rhot_data = None
        rhos_data = None

        ## clear intermediate files
        if settings['store_output_locally'] & settings['clear_output_zip_files']:
            for zfile in rhot_files+geom_files+rhos_files:
                if os.path.exists(zfile):
                    print('Removing {}'.format(zfile))
                    os.remove(zfile)

    return({gatts['acolite_type']:ofile})
