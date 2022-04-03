## def get gem from GEE
##
## written by Quinten Vanhellemont, RBINS
## 2021-02-26
## modifications: 2021-02-28 (QV) added to acolite, new default naming of gems
##                2022-03-07 (QV) added landsat collection keyword, added L9
##                2022-03-16 (QV) update for S2 PB4 offsets

def extract(st_lon, st_lat, sdate,
                           use_pid_name=False,
                           edate = None, # end date - if None sdate + 1 day will be used
                           st_name = None,
                           check_dates = None, # list of fyears to compare found images to
                           check_tiles = None,
                           max_diff_h = 3, # max difference between image time and fyear in check_dates
                           return_scene_list = False,
                           output = None,
                           gem_type = 'nc', # json
                           width = 100, override = False,
                           return_gem = False, return_iml = False,
                           landsat_collection = '02',
                           surface_reflectance = False,
                           sources = ['Landsat 5', 'Landsat 7','Landsat 8','Landsat 9', 'Sentinel-2'],
                           verbosity=0):

    import os, sys
    import json,bz2, dateutil.parser, datetime, re
    import numpy as np
    from pyproj import Proj

    import acolite as ac

    if check_tiles is not None:
        if type(check_tiles) is not list:
            check_tiles = list(check_tiles)

    ## set up ee
    try:
        import ee
        if False: ee.Authenticate()
        ee.Initialize()
    except:
        print(sys.exc_info()[0])
        print('Error initialising EarthEngine')
        return()

    ## make list of collections to search
    collections = []
    for source in sources:
        if source == 'Sentinel-2':
            if surface_reflectance:
                collections += ['COPERNICUS/S2_SR']
            else:
                collections += ['COPERNICUS/S2']
        if source == 'Landsat 8':
            if surface_reflectance:
                collections += ['LANDSAT/LC08/C{}/T1_L2'.format(landsat_collection), 'LANDSAT/LC08/C{}/T2_L2'.format(landsat_collection)]
            else:
                collections += ['LANDSAT/LC08/C{}/T1_TOA'.format(landsat_collection), 'LANDSAT/LC08/C{}/T2_TOA'.format(landsat_collection)]
        if source == 'Landsat 9':
            if surface_reflectance:
                collections += ['LANDSAT/LC09/C{}/T1_L2'.format(landsat_collection), 'LANDSAT/LC09/C{}/T2_L2'.format(landsat_collection)]
            else:
                collections += ['LANDSAT/LC09/C{}/T1_TOA'.format(landsat_collection), 'LANDSAT/LC09/C{}/T2_TOA'.format(landsat_collection)]
        if source == 'Landsat 5':
            if surface_reflectance:
                collections += ['LANDSAT/LT05/C{}/T1_L2'.format(landsat_collection), 'LANDSAT/LT05/C{}/T2_L2'.format(landsat_collection)]
            else:
                collections += ['LANDSAT/LT05/C{}/T1_TOA'.format(landsat_collection), 'LANDSAT/LT05/C{}/T2_TOA'.format(landsat_collection)]
        if source == 'Landsat 7':
            if surface_reflectance:
                collections += ['LANDSAT/LE07/C{}/T1_L2'.format(landsat_collection), 'LANDSAT/LE07/C{}/T2_L2'.format(landsat_collection)]
            else:
                collections += ['LANDSAT/LE07/C{}/T1_TOA'.format(landsat_collection), 'LANDSAT/LE07/C{}/T2_TOA'.format(landsat_collection)]

    ## check width
    width = min(width, 511) ## max allowed pixels is 512x512, since the box is centred on a pixel, the max we can ask is 511
    ## half box width
    hwidth = (width-1)/2

    ## set bounding dates
    ee_sdate=ee.Date(sdate)
    if edate is None: edate = (dateutil.parser.parse(sdate)+datetime.timedelta(days=1)).isoformat()[0:10]
    ee_edate=ee.Date(edate)

    ## make point geometry
    pt = ee.Geometry.Point([st_lon, st_lat])

    ## make location name
    loc_name = '_'.join(['{}'.format(abs(st_lat)).replace('.','N' if st_lat >=0 else 'S'),\
                     '{}'.format(abs(st_lon)).replace('.','E' if st_lon >=0 else 'W'), '{}px'.format(width)])
    if st_name is not None: loc_name = '{}_{}'.format(st_name, loc_name)
    if surface_reflectance: loc_name+='_SR'
    obase = '{}/{}'.format(output if output is not None else os.getcwd(), loc_name)

    ## search collections
    imColl = None
    for coll in collections:
        imC = ee.ImageCollection(coll).filterBounds(pt).filterDate(ee_sdate, ee_edate)
        if imColl is None:
            imColl = imC
        else:
            imColl = imColl.merge(imC)

    nimages = len(imColl.getInfo()['features'])
    if verbosity>0: print('Found {} image{}{}'.format(nimages, 's' if nimages !=1 else '', '' if st_name is None else ' for {}'.format(st_name)))
    if nimages == 0: return([])

    iml = imColl.toList(nimages).getInfo()
    if return_iml: return(iml)

    ## track scenes
    scene_list = []

    ## run extraction
    gemfiles = []
    iml = imColl.toList(nimages).getInfo()
    for ii in range(nimages):
        meta = iml[ii]['properties']

        if 'PRODUCT_ID' in meta: ## Sentinel-2 image
            fkey = 'PRODUCT_ID'
            pid = meta[fkey]
            satellite = pid[0:3]
            sensor = '{}_MSI'.format(satellite)
            bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']
            if surface_reflectance:
                bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B11', 'B12', 'SCL']
            scale_factor = 0.0001
            add_factor = 0

            ## for some PB4 products the add offset values are missing
            if meta['PROCESSING_BASELINE'][1]>='4': add_factor = -1000

            ## find band offsets introduces in PB4
            offsets = []
            for b in bands:
                if 'RADIO_ADD_OFFSET_{}'.format(b) in meta:
                    offsets.append(float(meta['RADIO_ADD_OFFSET_{}'.format(b)]))
            if len(offsets) > 0: add_factor = offsets[0]

            convert_counts = True ## S2 data still as counts on GEE
            target_res_default = 10
            tar_band = 'B2'
            dtime = [k for k in meta['GRANULE_ID'].split('_') if len(k) == 15][0]
            tile_name = '{}'.format(meta['MGRS_TILE'])
        elif 'LANDSAT_PRODUCT_ID' in meta: ## Landsat image
            fkey = 'LANDSAT_PRODUCT_ID'
            pid = meta[fkey]
            satellite = pid[0]+pid[3]
            if satellite == 'L5':
                sensor = '{}_TM'.format(satellite)
                bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7']
                if surface_reflectance:
                    bands = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']
            elif satellite == 'L7':
                sensor = '{}_ETM'.format(satellite)
                bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6_VCID_1', 'B6_VCID_2', 'B7', 'B8']
                if surface_reflectance:
                    bands = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7']
            elif satellite in ['L8', 'L9']:
                sensor = '{}_OLI'.format(satellite)
                bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11']
                if surface_reflectance:
                    bands = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10']

            if (landsat_collection == '02') & (not surface_reflectance): bands += ['SAA', 'SZA', 'VAA', 'VZA']

            refl_mult = np.unique([meta[k] for k in meta if 'REFLECTANCE_MULT' in k])[0]
            refl_add = np.unique([meta[k] for k in meta if 'REFLECTANCE_ADD' in k])[0]
            scale_factor = 1*refl_mult
            add_factor = 1*refl_add
            convert_counts = False ## Landsat data already converted on GEE

            if surface_reflectance:
                convert_counts = True
                scale_factor = 0.0000275
                add_factor = -0.2
                scale_factor_thermal = 0.00341802
                add_factor_thermal = -0.2149

            target_res_default = 30
            tar_band = bands[0]
            dtime = meta['DATE_ACQUIRED']+'T'+meta['SCENE_CENTER_TIME']

            row = '{}'.format(meta['WRS_ROW']).zfill(3)
            path = '{}'.format(meta['WRS_PATH']).zfill(3)
            tile_name = '{}{}'.format(path, row)

        if check_tiles is not None:
            if tile_name not in check_tiles: continue

        ## parse date time
        dt = dateutil.parser.parse(dtime)
        if use_pid_name:
            bname = '{}'.format(pid)
        else:
            bname = '{}_{}_{}'.format(sensor, dt.strftime('%Y_%m_%d_%H_%M_%S'), tile_name)

        ## preferably use NetCDF gems as they can be processed with acolite
        if gem_type == 'json':
            gemfile = '{}/{}_GEM.bz2'.format(obase, bname)
        elif gem_type == 'nc':
            gemfile = '{}/{}_GEM_L1R.nc'.format(obase, bname)
            if surface_reflectance:
                gemfile = '{}/{}_GEM_L2A.nc'.format(obase, bname)

        ## check time difference
        ## dates is a list of fractional years
        ## only keep scenes within max_diff_h of an element of dates
        if check_dates is not None:
            yf = ac.shared.isodate_to_yday(dt, return_yf=True)
            ylen = (datetime.datetime(dt.year+1, 1, 1) - datetime.datetime(dt.year, 1, 1)).days
            maxd = ((1/ylen) / 24) * max_diff_h
            mind = np.nanmin(np.abs(check_dates-yf))
            if mind>maxd:
                if verbosity > 1: print('Skipping {} {}>{}'.format(pid, mind, maxd))
                continue

        ## no need to transfer data if we just want list of scenes
        if return_scene_list:
            scene_list.append(meta[fkey])
            continue

        ## print some info on current image
        if verbosity>2:
            print(ii+1, nimages, datetime.datetime.now().isoformat()[0:19])
            print(pid)
            print(dtime)

        ## get data
        if not os.path.exists(gemfile) or (override):
            #dt = dateutil.parser.parse(dtime)
            isodate = dt.isoformat()
            fy = ac.shared.isodate_to_yday(isodate)
            fyear = fy[1]+fy[2]

            ## select product
            im = imColl.filter(ee.Filter.eq(fkey, pid)).first()
            properties = im.propertyNames().getInfo()

            ## get geometry
            if sensor in ['S2A_MSI', 'S2B_MSI']:
                saa = meta['MEAN_SOLAR_AZIMUTH_ANGLE']#im.get('MEAN_SOLAR_AZIMUTH_ANGLE').getInfo()
                sza = meta['MEAN_SOLAR_ZENITH_ANGLE']#im.get('MEAN_SOLAR_ZENITH_ANGLE').getInfo()
                geometry = {}
                for b in bands:
                    if 'QA' in b: continue
                    if 'SCL' in b: continue
                    if 'MEAN_INCIDENCE_AZIMUTH_ANGLE_{}'.format(b) in meta:
                        vaa = meta['MEAN_INCIDENCE_AZIMUTH_ANGLE_{}'.format(b)]#im.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_{}'.format(b)).getInfo()
                    else: continue
                    if 'MEAN_INCIDENCE_ZENITH_ANGLE_{}'.format(b) in meta:
                        vza = meta['MEAN_INCIDENCE_ZENITH_ANGLE_{}'.format(b)]#im.get('MEAN_INCIDENCE_ZENITH_ANGLE_{}'.format(b)).getInfo()
                    else: continue
                    geometry[b] = {'vza':vza, 'vaa':vaa}
                vza = np.nanmean([geometry[b]['vza'] for b in geometry])
                vaa = np.nanmean([geometry[b]['vaa'] for b in geometry])
            else:
                saa = meta['SUN_AZIMUTH']#im.get('SUN_AZIMUTH').getInfo()
                sza = 90 - meta['SUN_ELEVATION']#im.get('SUN_ELEVATION').getInfo()
                vaa = 0

                ## get scene bounding box
                x = [i[0] for i in meta['system:footprint']['coordinates']]
                y = [i[1] for i in meta['system:footprint']['coordinates']]

                ## extract corners
                wi = x.index(min(x)) ## west corner
                ei = x.index(max(x)) ## east corner
                ni = y.index(max(y)) ## north corner
                si = y.index(min(y)) ## south corner

                ## compute scene azimuth
                vaa = (ac.shared.azimuth_two_points(x[ni], y[ni], x[ei], y[ei]) +
                       ac.shared.azimuth_two_points(x[wi], y[wi], x[si], y[si]))/2
                vza = 0

            ## target band for subsetting projection
            proj = im.select(tar_band).projection().getInfo()
            if 'crs' in proj:
                p = ee.Projection(proj['crs'], proj['transform'])
            elif 'wkt' in proj:
                    p = ee.Projection(proj['wkt'], proj['transform'])
            scale = p.nominalScale().getInfo()

            ## get pixel coordinates in x/y
            tmp = ee.Image.clip(im.pixelCoordinates(im.select(tar_band).projection()),pt)
            ret = ee.Image.reduceRegion(tmp, ee.Reducer.toList()).getInfo()
            im_x, im_y = ret['x'][0], ret['y'][0]

            ## this part of the code extracts a single pixel
            if False:
                ## get pixel coordinates in lon/lat
                tmp = ee.Image.clip(im.pixelLonLat(),pt)
                ret = ee.Image.reduceRegion(tmp, ee.Reducer.toList(), scale=scale).getInfo()
                im_lon, im_lat = ret['longitude'][0], ret['latitude'][0]

                ## get pixel values
                pix = ee.Image.sampleRectangle(im, pt, defaultValue=0)
                pix_data = pix.getInfo()

                ## extract data from pix_data
                pdata = {b:np.asarray(pix_data['properties'][b]) for b in bands}
                if convert_counts:
                    for b in pdata:
                        if 'SCL' in b: continue
                        mask = np.where(pdata[b] == 0)
                        if sensor in ['S2A_MSI', 'S2B_MSI']:
                            pdata[b] = (pdata[b]+add_factor)*scale_factor
                        else:
                            if 'ST_' in b:
                                pdata[b] = (pdata[b]*scale_factor_thermal)+add_factor_thermal
                            else:
                                pdata[b] = (pdata[b]*scale_factor)+add_factor
                        pdata[b][mask] = np.nan

            ## use pixel coordinates from image to make new subset
            mins = ee.List([im_x-hwidth-0.5, im_y-hwidth-0.5])
            maxs = ee.List([im_x+hwidth+0.5, im_y+hwidth+0.5])
            rect = ee.Geometry.Rectangle(mins.cat(maxs), p, True, False).transform("EPSG:4326")

            ## reproject subset to target resolution
            ## so we have all bands at same dimensions
            rep = im.reproject(p, None, target_res_default * 1.0)
            box = ee.Image.sampleRectangle(rep, rect, defaultValue=0)
            box_data = box.getInfo()

            ## extract data from box_data
            data = {b:np.asarray(box_data['properties'][b]) for b in bands}
            for b in data:
                mask = np.where(data[b] == 0)
                if b in ['SAA', 'SZA', 'VAA', 'VZA']:
                    data[b] = data[b].astype(np.float32) / 100
                else:
                    if 'SCL' in b: continue
                    if convert_counts: ## convert from DN if needed
                        if sensor in ['S2A_MSI', 'S2B_MSI']:
                            data[b] = (data[b]+add_factor)*scale_factor
                        else:
                            if 'ST_' in b:
                                data[b] = (data[b]*scale_factor_thermal)+add_factor_thermal
                            else:
                                data[b] = (data[b]*scale_factor)+add_factor
                    data[b] = data[b].astype(np.float32)
                data[b][mask] = np.nan

            ## update average geometry
            if 'SAA' in data: saa = np.nanmean(data['SAA'])
            if 'VAA' in data: vaa = np.nanmean(data['VAA'])
            if 'SZA' in data: saa = np.nanmean(data['SZA'])
            if 'VZA' in data: vaa = np.nanmean(data['VZA'])
            if ('VAA' in data) & ('SAA' in data):
                data['RAA'] = np.abs(data['SAA'] - data['VAA'])
                data['RAA'][data['RAA']>=180] = np.abs(360-data['RAA'][data['RAA']>=180])
            raa = saa-vaa
            if raa < 0: raa = np.abs(raa)
            if raa >= 180: raa = np.abs(360 - raa)

            ## get lon and lat
            mp = im.pixelLonLat().reproject(p, None, target_res_default * 1.0)
            ret = ee.Image.sampleRectangle(mp, rect, defaultValue=0)
            ret_val = ret.getInfo()
            data['lon'] = np.asarray(ret_val['properties']['longitude'])
            data['lat'] = np.asarray(ret_val['properties']['latitude'])

            ## flip the data so north is up in the NetCDF
            flip = True
            if flip:
                for k in data: data[k] = np.flipud(data[k])

            ## store projection
            minx, miny = mins.getInfo()
            maxx, maxy = maxs.getInfo()

            ## figure out projection
            ## try parsing crs directly with proj
            ## for some versions of pyproj this fails, then get utm zone from epsg
            t = proj['transform']
            try:
                if 'crs' in proj:
                    p = Proj(proj['crs'])
                elif 'wkt' in proj:
                    p = Proj(proj['wkt'])
                proj4_string = p.crs.to_proj4()
            except:
                epsg = int(proj['crs'].split(':')[1])
                datum = 'WGS84'
                if 32600 < epsg <= 32660:
                    zone = epsg - 32600
                    proj4_list = ['+proj=utm',
                                  '+zone={}'.format(zone),
                                  '+datum={}'.format(datum),
                                  '+units=m',
                                  '+no_defs ']
                if 32700 < epsg <= 32760:
                    zone = epsg - 32700
                    proj4_list = ['+proj=utm',
                                  '+zone={}'.format(zone),
                                  '+south',
                                  '+datum={}'.format(datum),
                                  '+units=m',
                                  '+no_defs ']
                proj4_string = ' '.join(proj4_list)
                p = Proj(proj4_string)
                if verbosity>1:
                    print('Proj did not convert {}'.format(proj['crs']))
                    print('Assuming WGS84 with UTM')
                    print(proj4_string)

            ## compute projection info
            pixel_size = t[0], t[4]
            xdim = maxx - minx
            ydim = maxy - miny
            x0 = t[2] + (minx * pixel_size[0])
            x1 = t[2] + (maxx * pixel_size[0])
            y0 = t[5] + (miny * pixel_size[1])
            y1 = t[5] + (maxy * pixel_size[1])
            xrange = [x0, x1]
            yrange = [y0, y1]

            ## store projection info in dct to be stored in nc gatts
            is_utm = '+proj=utm' in proj4_string
            zone = -999
            m = re.search('\+zone=(.+?) ', proj4_string)
            if m: zone = int(m.group(1))
            dct = {'p': p,
                   #'epsg':  p.crs.to_epsg(),
                   'xrange': xrange, 'yrange': yrange,
                   'proj4_string':proj4_string,
                   'xdim':xdim, 'ydim':ydim,
                   'dimensions': (ydim, xdim),
                   'pixel_size': pixel_size, 'zone':zone}

            ## store matchup data
            gem = {'data':data, 'meta':meta, 'dct':dct,
                   'isodate':isodate, 'fyear': fyear,
                   'pid':pid, 'sensor':sensor,
                   'vza':vza, 'sza':sza, 'raa':raa, 'saa':saa, 'vaa': vaa}
            if surface_reflectance: gem['level'] = 'L2A'
            if return_gem: return(gem)

            if gem_type == 'json':
                gem_json_write(gemfile, gem)
            elif gem_type == 'nc':
                ac.gem.write(gemfile, gem)
        else:
            if verbosity>1: print('{} exists'.format(gemfile))

        gemfiles.append(gemfile)
    if return_scene_list:
        return(scene_list)
    else:
        return(gemfiles)
