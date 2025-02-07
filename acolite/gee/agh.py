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
##                2022-12-25 (QV) added SR option
##                2022-12-27 (QV) fixed SR computation (for ST data) and added hybrid/offline TACT run
##                2023-01-02 (QV) updated projection handling and scene centre lat/lon computation
##                2023-02-01 (QV) added extra parameters output for L2 ST data
##                2023-06-21 (QV) new version using computePixels, previous version renamed to agh_old
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-06-15 (QV) added to S2_SR_HARMONIZED
##                2025-02-05 (QV) added S2C_MSI

def agh(image, imColl, rsrd = {}, lutd = {}, luti = {}, settings = {}):
    import os, datetime, dateutil.parser, requests, json, time
    import acolite as ac
    from acolite import gee ## currently not imported in main acolite

    import numpy as np
    from pyproj import Proj
    from osgeo import gdal,osr
    gdal.UseExceptions()

    import ee
    #ee.Authenticate() ## assume ee use is authenticated in current environment
    #ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')

    uoz = settings['uoz_default']
    uwv = settings['uwv_default']
    pressure = settings['pressure']
    wind = settings['wind']

    rhop_par = settings['rhop_par']
    sel_par = settings['sel_par']

    reptran = settings['reptran']
    source = settings['source']
    emissivity = settings['emissivity']

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
    if ('LANDSAT' in fkey) & (pid[0:4] in ['LT08', 'LT09']): tar_band = 'B10' ## TIRS only data
    if ('LANDSAT' in fkey) & (pid[0:4] in ['LM01', 'LM02','LM03']): tar_band = 'B4'
    if ('LANDSAT' in fkey) & (settings['surface_reflectance']): tar_band = 'SR_B2' ## Landsat SR data

    ## check whether TACT is possible
    if settings['run_hybrid_tact']:
        if ('LANDSAT' not in fkey) :
            print('Not Landsat, no hybrid TACT processing possible')
            return()
        if (pid[0:4] in ['LO08', 'LO09']):
            print('OLI only data, no hybrid TACT processing possible')
            return()

    ## output file names for combined file
    ## used to determine NetCDF name, also for Google Drive outputs
    ## now .tif rather than .zip
    ext = ''
    if len(rname) >= 0: ext = '_{}'.format(rname)
    rhot_file = pid+'_rhot'+ext
    rhot_file_local = '{}/{}.tif'.format(output,rhot_file)
    rhos_file = pid+'_rhos'+ext
    rhos_file_local = '{}/{}.tif'.format(output,rhos_file)
    geom_file = pid+'_geom'+ext
    geom_file_local = '{}/{}.tif'.format(output,geom_file)
    sr_file = pid+'_sr'+ext
    sr_file_local = '{}/{}.tif'.format(output,sr_file)
    st_file = pid+'_st'+ext
    st_file_local = '{}/{}.tif'.format(output,st_file)
    sp_file = pid+'_sp'+ext
    sp_file_local = '{}/{}.tif'.format(output,sp_file)

    ## file type for file name
    file_type = 'L1R_GEE'
    if settings['run_hybrid_dsf']: file_type = 'L2R_GEE_HYBRID'
    if settings['run_hybrid_tact']:
        if settings['run_hybrid_dsf']:
            file_type = 'L2R_ST_GEE_HYBRID'
        else:
            file_type = 'ST_GEE_HYBRID'
    if settings['surface_reflectance']: file_type = 'SR_GEE'

    ## select product
    i = imColl.filter(ee.Filter.eq(fkey, pid)).first()

    ## get projection
    proj = i.select(tar_band).projection().getInfo()
    if 'crs' in proj:
        proj_crs = proj['crs']
    elif 'wkt' in proj:
        proj_crs = proj['wkt']
    p = ee.Projection(proj_crs, proj['transform'])
    scale = p.nominalScale().getInfo()
    if settings['output_scale'] is not None: scale = settings['output_scale']

    # #if (limit is not None) & (settings['subset_aot']) & (settings['run_hybrid_tact']):
    # if (limit is not None):
    #     if settings['strict_subset']:
    #         ## determin strict lat/lon rectangle box
    #         region = ee.Geometry.BBox(limit[1], limit[0], limit[3], limit[2])
    #     else:
    #         ## determine image bounding box
    #         imx = []
    #         imy = []
    #         for ii in [[1,0], [1,2], [3,0], [3,2]]:
    #             ## make point geometry
    #             pt = ee.Geometry.Point([limit[ii[0]], limit[ii[1]]])
    #
    #             ## get pixel coordinates in x/y
    #             tmp = ee.Image.clip(i.pixelCoordinates(i.select(tar_band).projection()),pt)
    #             ret = ee.Image.reduceRegion(tmp, ee.Reducer.toList()).getInfo()
    #             imx.append(ret['x'][0])
    #             imy.append(ret['y'][0])
    #
    #         ## use pixel coordinates from image to make new subset
    #         eesub = ee.List([min(imx), min(imy), max(imx), max(imy)])
    #         region = ee.Geometry.Rectangle(eesub, p, True, False)
    #
    # ## subset here if local aot is to be computed
    # #if (limit is not None) & (settings['subset_aot']) & (settings['run_hybrid_tact']): i = i.clip(region)
    # if (limit is not None) & (settings['subset_aot']): i = i.clip(region)

    ## get image info
    im = i.getInfo()
    bands = [p['id'] for p in im['bands'] if p['id'][0]]
    ## get projection info and region dimensions
    for b in im['bands']:
        if b['id'] == tar_band:
            odim = b['dimensions']
            if 'origin' in b:
                 origin = b['origin']
                 print('Region dimensions {}'.format(odim))
                 print('Region origin {}'.format(origin))
            else:
                 origin = 0, 0
                 print('Scene dimensions {}'.format(odim))
                 print('Scene origin {}'.format(origin))
            crs = b['crs']
            crs_transform = b['crs_transform']
            transform = b['crs_transform']
            print(crs)
            print(crs_transform)

    ## set up projection
    # prj = osr.SpatialReference()
    # prj.ImportFromEPSG(int(crs.split(':')[-1]))
    # Wkt = prj.ExportToWkt()
    # wp = Proj(Wkt)
    ## find scene center
    # nx = origin[0] + odim[0]/2
    # ny = origin[1] + odim[1]/2
    nx = origin[0] + odim[0]/2
    ny = origin[1] + odim[1]/2
    wp = Proj(proj_crs)
    mlon, mlat = wp(crs_transform[2]+nx*crs_transform[0], crs_transform[5]+ny*crs_transform[4], inverse=True)
    ll = {'longitude': mlon, 'latitude': mlat}
    print('Scene centre: {:.5f}E {:.5f}N'.format(ll['longitude'], ll['latitude']))

    if 'PRODUCT_ID' in im['properties']: ## Sentinel-2 image
        fkey = 'PRODUCT_ID'
        pid = im['properties'][fkey]
        dtime = [k for k in im['properties']['GRANULE_ID'].split('_') if len(k) == 15][0]
        tile_name = '{}'.format(im['properties']['MGRS_TILE'])
        satellite = pid[0:3]
        sensor = '{}_MSI'.format(satellite)
        satellite_sensor = '{}_MSI'.format(satellite)
        scale_factor = 0.0001
        add_factor = 0
        if (im['properties']['PROCESSING_BASELINE'][1]>='4') & ~(('S2_HARMONIZED' in im['id']) | ('S2_SR_HARMONIZED' in im['id'])) :
            add_factor = -1000

    elif 'LANDSAT_PRODUCT_ID' in im['properties']: ## Landsat image
        fkey = 'LANDSAT_PRODUCT_ID'
        pid = im['properties'][fkey]
        dtime = im['properties']['DATE_ACQUIRED']+'T'+im['properties']['SCENE_CENTER_TIME']
        row = '{}'.format(im['properties']['WRS_ROW']).zfill(3)
        path = '{}'.format(im['properties']['WRS_PATH']).zfill(3)
        tile_name = '{}{}'.format(path, row)
        satellite = pid[0]+pid[3]
        satellite_sensor = None

        if 'MSS' in im['properties']['SENSOR_ID']:
            sensor = '{}_MSS'.format(satellite)
        else:
            if satellite == 'L4':
                sensor = 'L5_TM' ## use same LUT as for L5
                satellite_sensor = '{}_TM'.format(satellite)
            elif satellite == 'L5':
                sensor = '{}_TM'.format(satellite)
            elif satellite == 'L7':
                sensor = '{}_ETM'.format(satellite)
            elif satellite in ['L8', 'L9']:
                sensor = '{}_OLI'.format(satellite)

        #print(im['properties'])

        if satellite_sensor is None: satellite_sensor = '{}'.format(sensor)
        scale_factor = 1
        add_factor = 0
        if settings['surface_reflectance']:
            scale_factor = 2.75e-05
            add_factor = -0.2
            thermal_scale_factor = 0.00341802
            thermal_add_factor = 149

    ## parse datetime
    dt = dateutil.parser.parse(dtime)

    if settings['run_hybrid_tact']:
        max_date = (datetime.datetime.now() - datetime.timedelta(days=91)).isoformat()
        if dt.isoformat() > max_date:
            print('File too recent for TACT with {} profiles: after {}'.format('era5', max_date))
            #print('Run with tact_profile_source=gdas1 or tact_profile_source=ncep.reanalysis2 for NRT processing')
            return()

    if settings['use_scene_name']:
        ofile = rhot_file_local.replace('_rhot', '_{}'.format(file_type)).replace('.tif', '.nc')
    else:
        obase  = '{}_{}_{}_{}{}'.format(satellite_sensor,  dt.strftime('%Y_%m_%d_%H_%M_%S'), tile_name, file_type, ext)
        ofile = '{}/{}.nc'.format(os.path.dirname(rhot_file_local), obase)

    if os.path.exists(ofile) & (settings['override'] is False):
        print('AGH file {} exists, set override=True to replace'.format(ofile))
        return()

    if sensor not in rsrd:
        ## get sensor defaults
        setd = ac.acolite.settings.parse(sensor)
        if 'rsr_version' in setd:
            lut_sensor = '{}_{}'.format(sensor, setd['rsr_version'])
        else:
            lut_sensor = '{}'.format(sensor)

        print('Loading RSRs: {}'.format(sensor))
        rsrd[sensor] = ac.shared.rsr_dict(sensor=lut_sensor)[lut_sensor]

    ## select product again
    ## cropped product seems to give empty tiles
    ##i = imColl.filter(ee.Filter.eq(fkey, pid)).first()

    ## if processing rhot
    if not settings['surface_reflectance']:
        ## get central lon and lat
        #if settings['ancillary_data'] | settings['run_hybrid_tact']:
        #    mp = i.pixelLonLat()
        #    ll = mp.reduceRegion(crs=proj_crs, crsTransform=proj['transform'], geometry=(i.geometry() if limit is None else region), \
        #                         reducer= ee.Reducer.mean(), bestEffort=True).getInfo()

        ## get ancillary
        if settings['ancillary_data']:
            print('Getting ancillary for scene centre location.')
            # anc = ac.ac.ancillary.get(dt, ll['longitude'], ll['latitude'])
            #
            # ## overwrite the defaults
            # if ('ozone' in anc): uoz = anc['ozone']['interp']/1000. ## convert from MET data
            # if ('p_water' in anc): uwv = anc['p_water']['interp']/10. ## convert from MET data
            # if ('z_wind' in anc) & ('m_wind' in anc) & (wind is None):
            #     wind = ((anc['z_wind']['interp']**2) + (anc['m_wind']['interp']**2))**0.5
            # if ('press' in anc): pressure = anc['press']['interp']
            # print(uoz, uwv, wind, pressure)

            anc = ac.ac.ancillary.get(dt, ll['longitude'], ll['latitude'])
            if ('uoz' in anc): uoz = anc['uoz']
            if ('uwv' in anc): uwv = anc['uwv']
            if ('wind' in anc): wind = anc['wind']
            if ('pressure' in anc): pressure = anc['pressure']
            print(uoz, uwv, wind, pressure)

        if settings['run_hybrid_tact']:
            emissivity_file = '{}/{}/emissivity_{}.json'.format(ac.config['data_dir'], 'TACT', emissivity)
            #print(emissivity_file)
            print('Getting TACT parameters for scene centre location.')
            thermal_sensors = {'5':'L5_TM', '7':'L7_ETM', '8':'L8_TIRS', '9':'L9_TIRS'}
            thermal_bands = {'5':'6', '7':['6_VCID_1','6_VCID_2'], '8':'10', '9':'10'}
            thermal_sensor = thermal_sensors[sensor[1]]
            #thermal_band = thermal_bands[sensor[1]]

            ## radiative transfer
            thd, simst, lonc, latc = ac.tact.tact_limit(dtime,
                                                        lon = ll['longitude'],
                                                        lat = ll['latitude'],
                                                        satsen = thermal_sensor,
                                                        reptran = reptran, source = source)
            #print(thd)
            em = json.load(open(emissivity_file, 'r'))
            #bem = em[thermal_sensor][thermal_band]
            #btau = thd['tau{}'.format(thermal_band)]
            #bLu = thd['Lu{}'.format(thermal_band)]
            #bLd = thd['Ld{}'.format(thermal_band)]

        ## load luts
        if (sensor not in lutd) & (settings['run_hybrid_dsf']):
            print('Loading atmosphere LUTs: {}'.format(sensor))
            lutd[sensor] = ac.aerlut.import_luts(sensor=sensor)

        if (sensor not in luti) & (settings['run_hybrid_dsf']) & (settings['glint_correction']):
            print('Loading interface LUTs: {}'.format(sensor))
            luti[sensor] = ac.aerlut.import_rsky_luts(models=[1,2], lutbase='ACOLITE-RSKY-202102-82W', sensor=sensor)

    ## get average geometry
    ## will be updated later if SAA/SZA/VAA/VZA data are available
    geometry = {}
    if sensor in ['S2A_MSI', 'S2B_MSI', 'S2C_MSI']:
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

    ## get either SR or rhot
    if settings['surface_reflectance']:
        info = i.getInfo()
        #print(info)
        ## make sr dataset
        i_sr = None
        for b in rsrd[sensor]['rsr_bands']:
            bname = 'B{}'.format(b)
            if sensor[0] == 'L': bname='SR_{}'.format(bname)
            try:
                if sensor[0] == 'L':
                    sr = i.expression('(DN * {}) + {}'.format(scale_factor,add_factor), {'DN': i.select(bname)})
                else:
                    sr = i.expression('(DN + {}) * {}'.format(add_factor,scale_factor), {'DN': i.select(bname)})
                tmp = sr.getInfo()
            except:
                continue
            print('sr {}'.format(bname))
            if i_sr is None:
                i_sr = ee.Image(sr)
            else:
                i_sr = i_sr.addBands(sr)

        ## add ST in thermal bands
        if sensor[0] == 'L':
            if sensor[1] in ['5', '7']:
                bnames = ['ST_B6']
            if sensor[1] in ['8', '9']:
                bnames = ['ST_B10'] # 'ST_B11' B11 data is not produced?
            for bname in bnames:
                print(sensor, bname)
                try:
                    st = i.expression('(DN * {}) + {} '.format(thermal_scale_factor,thermal_add_factor), {'DN': i.select(bname)})
                    tmp = st.getInfo()
                except:
                    continue
                i_sr = i_sr.addBands(st)

        ## SR band info
        sr_info = i_sr.getInfo()
        obands_sr = [ib['id'] for ib in sr_info['bands']]
        print(obands_sr)

        ## add extra par for ST
        st_par = None
        st_par_bands = []
        if 'PROCESSING_LEVEL' in im['properties']:
            if (sensor[0] == 'L') & (settings['store_sp']) & (im['properties']['PROCESSING_LEVEL'] == 'L2SP'):
                sp_scale = {'ST_B6': 0.00341802, 'ST_B10': 0.00341802,
                            'ST_ATRAN': 0.0001, 'ST_CDIST': 0.01,
                            'ST_DRAD': 0.001, 'ST_EMIS': 0.0001,
                            'ST_EMSD': 0.0001, 'ST_QA': 0.01,
                            'ST_TRAD': 0.001,'ST_URAD': 0.001}

                sp_offset = {'ST_B6': 149, 'ST_B10': 149,
                             'ST_ATRAN': 0, 'ST_CDIST': 0,
                             'ST_DRAD': 0, 'ST_EMIS': 0,
                             'ST_EMSD': 0,'ST_QA': 0,
                             'ST_TRAD': 0,'ST_URAD': 0}

                for b in im['bands']:
                    bname = b['id']
                    if bname[0:2] == 'ST':
                        print(bname)
                        sp = i.expression('(DN * {}) + {} '.format(sp_scale[bname],sp_offset[bname]),
                                          {'DN': i.select(bname)})
                        if st_par is None:
                            st_par = ee.Image(sp)
                        else:
                            st_par = st_par.addBands(sp)
                        st_par_bands.append(bname)
                    else:
                        continue
    else:
        ## make rhot dataset
        i_rhot = None
        for b in rsrd[sensor]['rsr_bands']:
            bname = 'B{}'.format(b)
            print('rhot {}'.format(bname))
            if 'MSS' in sensor:
                add_factor = im['properties']['REFLECTANCE_ADD_BAND_{}'.format(b)]
                scale_factor = im['properties']['REFLECTANCE_MULT_BAND_{}'.format(b)]
                if i_rhot is None: i = i.toFloat()

            ## gas and path corrected rhot
            rhot = i.expression('(DN + {}) * {}'.format(add_factor, scale_factor), {'DN': i.select(bname)})
            if 'MSS' in sensor:
                mus = np.cos(np.radians(90-im['properties']['SUN_ELEVATION']))
                rhot = rhot.expression('rhot/{}'.format(mus), {'rhot': rhot.select(bname)})

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
        print(obands_rhot)

        print('Getting geometry and gas transmittance')
        ## get geometry percentiles
        geom_percentile = 50
        prc = i.reduceRegion(reducer= ee.Reducer.percentile([geom_percentile]), bestEffort=True, maxPixels=1e13).getInfo() #
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
            prc = i_rhot.reduceRegion(reducer= ee.Reducer.percentile(percentiles), bestEffort=True, maxPixels=1e13).getInfo()
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
            print('Performing DSF atmospheric correction')
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

        ## run TACT on this data
        if settings['run_hybrid_tact']:
            ## do atmospheric correction
            print('Performing TACT atmospheric correction')
            i_st = None
            for bname in thermal_bands:
                b = bname.replace('B', '')

                ## K constants
                k1 = float(im['properties']['K1_CONSTANT_BAND_{}'.format(b)])
                k2 = float(im['properties']['K2_CONSTANT_BAND_{}'.format(b)])
                if thermal_sensor == 'L7_ETM': b = '6' ## avoid VCID distinction

                print(bname, k1, k2)
                bem = em[thermal_sensor][b]
                btau = thd['tau{}'.format(b)]
                bLu = thd['Lu{}'.format(b)]
                bLd = thd['Ld{}'.format(b)]

                ## Lt
                lt = i_rhot.expression('{k1}/(exp({k2}/BT)-1)'.format(k1=k1, k2=k2), \
                                        {'BT': i_rhot.select(bname)}).rename(bname)

                ## Ls
                ls = lt.expression('(((LT-{Lu})/{tau}) - ((1-{em})*{Ld}))/{em}'.format(em=bem, Lu=bLu, Ld=bLd, tau=btau), \
                                        {'LT': lt.select(bname)}).rename(bname)

                ## ST
                st = ls.expression('{k2}/(log(({k1}/LS)+1))'.format(k2=k2, k1=k1), \
                                        {'LS': ls.select(bname)}).rename(bname)

                if i_st is None:
                    i_st = ee.Image(st)
                else:
                    i_st = i_st.addBands(st)

            ## get band info for outputs
            st_info = i_st.getInfo()
            obands_st = [ib['id'] for ib in st_info['bands']]


    # ## check size and tile if necessary
    # tiled_transfer = False
    # tile_size = [606,606]
    # tile_size = [306,306]
    # tiles_grid = 1,1
    # if ((odim[0]*odim[1]*13)*9 > 50331648) or (odim[0]>tile_size[0]) or (odim[1]>tile_size[1]): #((odim[0]*odim[1]) > (tile_size[0]*tile_size[1])):
    #     print('Warning large dataset: {}x{}'.format(odim[0],odim[1]))
    #     tiled_transfer = True
    #     tiles_grid = int(np.ceil(odim[0]/tile_size[0])), int(np.ceil(odim[1]/tile_size[1]))
    # num_tiles = tiles_grid[0]*tiles_grid[1]
    #
    # ## identify tiles
    # if tiled_transfer:
    #     tiles = []
    #     tn = 0
    #     for ti in range(tiles_grid[0]):
    #         for tj in range(tiles_grid[1]):
    #             ti0 = ti*tile_size[0]
    #             ti1 = min((ti+1)*tile_size[0], odim[0])
    #             tj0 = tj*tile_size[1]
    #             tj1 = min((tj+1)*tile_size[1], odim[1])
    #             tiles.append([ti0, ti1, tj0, tj1, 't{}'.format(str(tn).zfill(3))])
    #             print('tile', tiles[-1])
    #             tn+=1
    # else:
    #     tiles = [[0,odim[0], 0,odim[1],'']]

    #tile_size = [606,606]
    tile_size = settings['tile_size']
    tiled_transfer = True

    ## get original image subset
    ## for getPixels computePixels grid
    ## can replace some subsetting code above as well, to be done
    if limit is None:
        porigin = 0.,0.
    else:
        ## project limit to image CRS
        bbox = np.asarray(wp((limit[1], limit[3], limit[3], limit[1]),
                            (limit[0], limit[0], limit[2], limit[2]), ))
        ## move to nearest pixel - 60 m works for both S2 and Landsat
        bbox -= bbox % 60
        ## ROI x and y ranges
        xrange = np.asarray((min(bbox[0]), max(bbox[0])))
        yrange = np.asarray((max(bbox[1]), min(bbox[1])))
        xdim = (xrange[1]-xrange[0])/transform[0]
        ydim = (yrange[1]-yrange[0])/transform[4]
        ## ROI x and y pixel ranges
        xprange = (xrange-transform[2]) / transform[0]
        yprange = (yrange-transform[5]) / transform[4]
        rdim = int(xprange[1]-xprange[0]), int(yprange[1]-yprange[0])
        ## origin of subset in image
        porigin = xprange[0], yprange[0]
    ## tiles in original pixel grid
    print(porigin)
    print(rdim)
    ptiles = []
    tn = 0
    tiles_grid = int(np.ceil(rdim[0]/tile_size[0])), int(np.ceil(rdim[1]/tile_size[1]))
    num_tiles = tiles_grid[0]*tiles_grid[1]
    for ti in range(tiles_grid[0]):
        for tj in range(tiles_grid[1]):
            ti0 = max(porigin[0], porigin[0]+ti*tile_size[0])
            ti1 = min(porigin[0]+(ti+1)*tile_size[0], porigin[0]+rdim[0])#, odim[0])
            tj0 = max(porigin[1], porigin[1]+tj*tile_size[1])
            tj1 = min(porigin[1]+(tj+1)*tile_size[1], porigin[1]+rdim[1])#, odim[1])
            ptiles.append([ti0, ti1, tj0, tj1, 't{}'.format(str(tn).zfill(3))])
            print('tile', ptiles[-1])
            tn+=1
    print(transform)
    print('Running download with {} tiles'.format(len(ptiles)))

    ## output data to drive
    if settings['store_output_google_drive']:
        ## set output config for rhot
        output_config = {'description': rhot_file, 'folder':'ACOLITE', 'scale': scale}

        ## test if setting crs works
        output_config['crs'] = proj_crs
        output_config['crs_transform'] = proj['transform']

        ## set output dir
        if settings['drive_output'] is not None: output_config['folder'] = settings['drive_output']
        if limit is not None:
            output_config['region'] = region
        else:
            output_config['region'] = i.geometry()

        if settings['surface_reflectance']:
            if settings['store_sr']:
                output_config['description'] = sr_file
                task_sr = ee.batch.Export.image.toDrive(i_sr.select(obands_sr), **output_config)
                print('Exporting to Drive {}'.format(output_config['description']))
                task_sr.start()
                status = gee.check_task(task_sr, task_sleep = settings['task_check_sleep'])
        else:
            if settings['store_rhot']:
                output_config['description'] = rhot_file
                task_rhot = ee.batch.Export.image.toDrive(i_rhot.select(obands_rhot), **output_config)
                print('Exporting to Drive {}'.format(output_config['description']))
                task_rhot.start()
                status = gee.check_task(task_rhot, task_sleep = settings['task_check_sleep'])

            if settings['store_rhos'] & settings['run_hybrid_dsf']:
                output_config['description'] = rhos_file
                task_rhos = ee.batch.Export.image.toDrive(i_rhos.select(obands_rhos), **output_config)
                print('Exporting to Drive {}'.format(output_config['description']))
                task_rhos.start()
                status = gee.check_task(task_rhos, task_sleep = settings['task_check_sleep'])

            if settings['store_st'] & settings['run_hybrid_tact']:
                output_config['description'] = st_file
                task_st = ee.batch.Export.image.toDrive(i_st.select(obands_st), **output_config)
                print('Exporting to Drive {}'.format(output_config['description']))
                task_st.start()
                status = gee.check_task(task_st, task_sleep = settings['task_check_sleep'])

        if settings['store_geom'] & (i_geom is not None):
            output_config['description'] = geom_file
            task_geom = ee.batch.Export.image.toDrive(i_geom.select(obands_geom), **output_config)
            print('Exporting to Drive {}'.format(output_config['description']))
            task_geom.start()
            status = gee.check_task(task_geom, task_sleep = settings['task_check_sleep'])

    ## end output to drive

    ## store data locally, with tiling if needed
    if settings['store_output_locally']:
        ## empty list to track local files
        rhot_files = []
        geom_files = []
        rhos_files = []
        sr_files = []
        st_files = []
        sp_files = []

        for ti, tile in enumerate(ptiles):
            ## output file names
            ext = ''
            if len(rname) >= 0: ext = '_{}'.format(rname)
            if tiled_transfer:
                tile_id = tile[4]
                ext+='_'+tile_id

            print(tile)
            tdim = tile[1]-tile[0], tile[3]-tile[2]
            print(tdim)
            tilex = transform[2] + tile[0] * transform[0]
            tiley = transform[5] + tile[2] * transform[4]
            #print(tile)
            print(tilex, tiley)

            ## set output names for this tile
            rhot_file_tile = pid+'_rhot'+ext
            rhot_file_tile_local = '{}/{}.tif'.format(output,rhot_file_tile)
            rhos_file_tile = pid+'_rhos'+ext
            rhos_file_tile_local = '{}/{}.tif'.format(output,rhos_file_tile)
            geom_file_tile = pid+'_geom'+ext
            geom_file_tile_local = '{}/{}.tif'.format(output,geom_file_tile)
            sr_file_tile = pid+'_sr'+ext
            sr_file_tile_local = '{}/{}.tif'.format(output,sr_file_tile)
            st_file_tile = pid+'_st'+ext
            st_file_tile_local = '{}/{}.tif'.format(output,st_file_tile)
            sp_file_tile = pid+'_sp'+ext
            sp_file_tile_local = '{}/{}.tif'.format(output,sp_file_tile)

            ## set up grid
            grid = {
                'dimensions': {
                    'width': tdim[0],
                    'height': tdim[1]
                },
                'affineTransform': {
                    'scaleX': transform[0],
                    'shearX': 0,
                    'translateX': tilex,
                    'shearY': 0,
                    'scaleY': transform[4],
                    'translateY': tiley,
                },
                'crsCode': proj_crs,
            }

            print(grid)

            ## store data
            ## write geometry (tile)
            if settings['store_geom'] & (i_geom is not None):
                ## set up request
                request = {'expression': (i_geom),'fileFormat': 'GEO_TIFF',
                            'grid': grid,'bandIds': obands_geom}
                ## get pixels
                tmp = ee.data.computePixels(request)
                print('geom', len(tmp))
                ## write data
                with open(geom_file_tile_local, 'wb') as f:
                    f.write(tmp)
                geom_files.append(geom_file_tile_local)

            ## write surface reflectance
            if settings['surface_reflectance']:
                if settings['store_sr']:
                    ## set up request
                    request = {'expression': (i_sr),'fileFormat': 'GEO_TIFF',
                                'grid': grid,'bandIds': obands_sr}
                    ## get pixels
                    tmp = ee.data.computePixels(request)
                    print('sr', len(tmp))
                    ## write data
                    with open(sr_file_tile_local, 'wb') as f:
                        f.write(tmp)
                    sr_files.append(sr_file_tile_local)
                if settings['store_sp'] & (st_par != None):
                    ## set up request
                    request = {'expression': (st_par),'fileFormat': 'GEO_TIFF',
                                'grid': grid,'bandIds': obands_sp}
                    ## get pixels
                    tmp = ee.data.computePixels(request)
                    print('sp', len(tmp))
                    ## write data
                    with open(sp_file_tile_local, 'wb') as f:
                        f.write(tmp)
                    sp_files.append(sp_file_tile_local)
            ## write rhot and rhos
            else:
                ## write rhot
                if settings['store_rhot']:
                    ## set up request
                    request = {'expression': (i_rhot),'fileFormat': 'GEO_TIFF',
                                'grid': grid,'bandIds': obands_rhot, 'workloadTag': tile[-1]}
                    ## get pixels
                    tmp = ee.data.computePixels(request)
                    print('rhot', len(tmp))
                    ## write data
                    with open(rhot_file_tile_local, 'wb') as f:
                        f.write(tmp)
                    rhot_files.append(rhot_file_tile_local)
                ## write rhos
                if settings['store_rhos'] & settings['run_hybrid_dsf']:
                    ## set up request
                    request = {'expression': (i_rhos),'fileFormat': 'GEO_TIFF',
                                'grid': grid,'bandIds': obands_rhos}
                    ## get pixels
                    tmp = ee.data.computePixels(request)
                    print('rhos', len(tmp))
                    ## write data
                    with open(rhos_file_tile_local, 'wb') as f:
                        f.write(tmp)
                    rhos_files.append(rhos_file_tile_local)
                ## write st
                if settings['store_st'] & settings['run_hybrid_tact']:
                    ## set up request
                    request = {'expression': (i_st),'fileFormat': 'GEO_TIFF',
                                'grid': grid,'bandIds': obands_st}
                    ## get pixels
                    tmp = ee.data.computePixels(request)
                    print('st', len(tmp))
                    ## write data
                    with open(st_file_tile_local, 'wb') as f:
                        f.write(tmp)
                    st_files.append(st_file_tile_local)
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
        sr_data = None
        st_data = None
        sp_data = None

        for ti, ptile in enumerate(ptiles):
            ## tile dimensions
            xp = ptile[1]-ptile[0]
            yp = ptile[3]-ptile[2]
            ## output tile position
            tile = int(ptile[0] - porigin[0]), int(ptile[1] - porigin[0]), \
                   int(ptile[2] - porigin[1]), int(ptile[3] - porigin[1])
            file_proj = None

            if settings['surface_reflectance']:
                sr_image_file = sr_files[ti]
                ## read sr
                if os.path.exists(sr_image_file):
                    file_proj = '{}'.format(sr_image_file)
                    ds = gdal.Open(sr_image_file)
                    if num_tiles == 1:
                        sr_data = ds.ReadAsArray()
                    else:
                        if sr_data is None: sr_data = np.zeros((ds.RasterCount, rdim[1], rdim[0]))+np.nan
                        sr_data[:, tile[2]:tile[3], tile[0]:tile[1]] = ds.ReadAsArray()
                    ds = None
                    sr_data[sr_data==0.0] = np.nan

                if settings['store_sp'] & (st_par != None):
                    sp_image_file = sp_files[ti]
                    ## read sp
                    if os.path.exists(sp_image_file):
                        file_proj = '{}'.format(sp_image_file)
                        ds = gdal.Open(sp_image_file)
                        if num_tiles == 1:
                            sp_data = ds.ReadAsArray()
                        else:
                            if sp_data is None: sp_data = np.zeros((ds.RasterCount, rdim[1], rdim[0]))+np.nan
                            sp_data[:, tile[2]:tile[3], tile[0]:tile[1]] = ds.ReadAsArray()
                        ds = None
                        #sr_data[sr_data==0.0] = np.nan
            else:
                ## read rhot data
                if len(rhot_files) == num_tiles:
                    rhot_image_file = rhot_files[ti]
                    ## read rhot
                    if os.path.exists(rhot_image_file):
                        file_proj = '{}'.format(rhot_image_file)
                        ds = gdal.Open(rhot_image_file)
                        if num_tiles == 1:
                            rhot_data = ds.ReadAsArray()
                        else:
                            if rhot_data is None: rhot_data = np.zeros((ds.RasterCount, rdim[1], rdim[0]))+np.nan
                            rhot_data[:, tile[2]:tile[3], tile[0]:tile[1]] = ds.ReadAsArray()
                        ds = None
                        rhot_data[rhot_data==0.0] = np.nan

                ## read rhos data
                if len(rhos_files) == num_tiles:
                    rhos_image_file = rhos_files[ti]
                    ## read rhos
                    if os.path.exists(rhos_image_file):
                        file_proj = '{}'.format(rhos_image_file)
                        ds = gdal.Open(rhos_image_file)
                        if num_tiles == 1:
                            rhos_data = ds.ReadAsArray()
                        else:
                            if rhos_data is None: rhos_data = np.zeros((ds.RasterCount, rdim[1], rdim[0]))+np.nan
                            rhos_data[:, tile[2]:tile[3],tile[0]:tile[1]] = ds.ReadAsArray()
                        ds = None
                    rhos_data[rhos_data==0.0] = np.nan

            ## read geom data
            if len(geom_files) == num_tiles:
                geom_image_file = geom_files[ti]
                ## read geom
                if os.path.exists(geom_image_file):
                    file_proj = '{}'.format(geom_image_file)
                    ds = gdal.Open(geom_image_file)
                    if num_tiles == 1:
                        geom_data = ds.ReadAsArray()
                    else:
                        if geom_data is None: geom_data = np.zeros((ds.RasterCount, rdim[1], rdim[0]))+np.nan
                        geom_data[:, tile[2]:tile[3],tile[0]:tile[1]] = ds.ReadAsArray()
                    ds = None

            ## read st data
            if len(st_files) == num_tiles:
                st_image_file = st_files[ti]
                print(st_image_file)
                ## read geom
                if os.path.exists(st_image_file):
                    file_proj = '{}'.format(st_image_file)
                    ds = gdal.Open(st_image_file)
                    if num_tiles == 1:
                        st_data = ds.ReadAsArray()
                    else:
                        if stdata is None: st_data = np.zeros((ds.RasterCount, rdim[1], rdim[0]))+np.nan
                        st_data[:, tile[2]:tile[3],tile[0]:tile[1]] = ds.ReadAsArray()
                    ds = None

            ## read file projection
            if ti == 0:
                dct = ac.shared.projection_read(file_proj)
            else:
                ## update dct
                dct_ = ac.shared.projection_read(file_proj)
                dct['yrange'] = max(dct['yrange'][0], dct_['yrange'][0]), min(dct['yrange'][1], dct_['yrange'][1])
                dct['xrange'] = min(dct['xrange'][0], dct_['xrange'][0]), max(dct['xrange'][1], dct_['xrange'][1])
                dct['ydim'] = int((dct['yrange'][1]-dct['yrange'][0])/dct['pixel_size'][1])
                dct['xdim'] = int((dct['xrange'][1]-dct['xrange'][0])/dct['pixel_size'][0])
                dct['dimensions'] = (dct['xdim'], dct['ydim'])

        ## image file in zip
        #rhot_image_file = '/vsizip/{}/{}.tif'.format(rhot_file_local, rhot_file)
        #rhos_image_file = '/vsizip/{}/{}.tif'.format(rhos_file_local, rhos_file)
        #geom_image_file = '/vsizip/{}/{}.tif'.format(geom_file_local, geom_file)

        ## gatts for output NetCDF
        gatts = {}
        gatts['acolite_type'] = 'l1r_gee'
        gatts['pid'] = pid
        gatts['sensor'] = sensor
        gatts['satellite_sensor'] = satellite_sensor
        gatts['isodate'] = dt.isoformat()
        ## add image metadata from properties
        property_keys = list(im['properties'].keys())
        property_keys.sort()
        for k in property_keys: gatts[k] = im['properties'][k]
        ## add geometry average
        for k in geometry: gatts[k] = geometry[k]
        #file_type = 'L1R_GEE'

        if settings['surface_reflectance']:
            gatts['acolite_type'] = 'sr_gee'
        else:
            if settings['run_hybrid_dsf']:
                gatts['acolite_type'] = 'l2r_gee_hybrid'
                gatts['pressure'] = pressure
                gatts['uwv'] = uwv
                gatts['uoz'] = uoz
                gatts['ac_aot_550'] = sel_aot
                gatts['ac_model'] = sel_lut
                gatts['ac_var'] = sel_val
            if settings['run_hybrid_tact']:
                if settings['run_hybrid_dsf']:
                    gatts['acolite_type'] = 'l2r_st_gee_hybrid'
                else:
                    gatts['acolite_type'] = 'st_gee_hybrid'

        ## output file names
        ext = ''
        if len(rname) >= 0: ext = '_{}'.format(rname)
        if settings['use_scene_name']:
            ofile = rhot_file_local.replace('_rhot', '_{}'.format(file_type)).replace('.tif', '.nc')
        else:
            obase  = '{}_{}_{}_{}{}'.format(gatts['satellite_sensor'],  dt.strftime('%Y_%m_%d_%H_%M_%S'), tile_name, file_type, ext)
            ofile = '{}/{}.nc'.format(os.path.dirname(rhot_file_local), obase)
        print(ofile)

        ## convert dct info into nc_projection and lat, lon
        #dct = ac.shared.projection_read(rhot_image_file)
        nc_projection = ac.shared.projection_netcdf(dct, add_half_pixel=True)
        lon, lat = ac.shared.projection_geo(dct, add_half_pixel=True)

        ## write lon
        ac.output.nc_write(ofile, 'lon', lon, attributes=gatts, new=new, double=True, nc_projection=nc_projection)
        if verbosity > 1: print('Wrote lon')
        lon = None

        ## write lat
        ac.output.nc_write(ofile, 'lat', lat, double=True)
        if verbosity > 1: print('Wrote lat')
        lat = None
        new = False

        if geom_data is not None:
            for bi, b in enumerate(['SAA', 'SZA', 'VAA', 'VZA']):
                ac.output.nc_write(ofile, b.lower(), geom_data[bi,:,:])
                if verbosity > 1: print('Wrote {}'.format(b.lower()))
            ## compute raa
            raa = np.abs(geom_data[0,:,:]-geom_data[2,:,:])
            raa[raa>180] = np.abs(360 - raa[raa>180])
            ac.output.nc_write(ofile, 'raa', raa)
            if verbosity > 1: print('Wrote raa')
            raa = None
            geom_data = None

        bii = 0 ## band index for SR files, which may not have all bands
        for bi, b in enumerate(rsrd[sensor]['rsr_bands']):
            wave_name = rsrd[sensor]['wave_name'][b]
            wave_nm = rsrd[sensor]['wave_nm'][b]
            att = {'band': b, 'wave_name': wave_name, 'wave_nm': wave_nm}

            if settings['surface_reflectance']:
                if sensor[0] == 'L':
                    if 'SR_B{}'.format(b) not in obands_sr: continue
                else:
                    if 'B{}'.format(b) not in obands_sr: continue
                sr_ds = 'rhos_sr_{}'.format(wave_name)
                ## write sr data
                if sr_data is not None:
                    ac.output.nc_write(ofile, sr_ds, sr_data[bii, :, :], dataset_attributes=att)
                    print('Wrote {}'.format(sr_ds))
                    bii += 1
            else:
                for k in ttg: att[k] = ttg[k][b]

                rhot_ds = 'rhot_{}'.format(wave_name)
                rhos_ds = 'rhos_{}'.format(wave_name)

                ## write rhot data
                if rhot_data is not None:
                    ac.output.nc_write(ofile, rhot_ds, rhot_data[bi, :, :], dataset_attributes=att)
                    print('Wrote {}'.format(rhot_ds))

                ## write rhos data
                if (att['tt_gas'] > min_tgas_rho) & (rhos_data is not None):
                    ac.output.nc_write(ofile, rhos_ds, rhos_data[bi, :, :], dataset_attributes=att)
                    print('Wrote {}'.format(rhos_ds))

        ## write Landsat L2 ST
        if settings['surface_reflectance'] & (sensor[0] == 'L'):
            for b in bnames:
                if b not in obands_sr: continue
                sr_ds = 'st{}'.format(b.replace('ST_B', ''))
                att = {'band': b}
                ## write st data
                if sr_data is not None:
                    ## mask out of scene data
                    cur_data = sr_data[bii, :, :] * 1.0
                    cur_data[cur_data<=thermal_add_factor] = np.nan
                    ac.output.nc_write(ofile, sr_ds, cur_data, dataset_attributes=att)
                    print('Wrote {} ({})'.format(sr_ds, cur_data.shape))
                    cur_data = None
                    bii += 1

            ## output extra st parameters
            if len(st_par_bands) != 0:
                for bi, b in enumerate(st_par_bands):
                    att = {'band': b}
                    if sp_data is not None:
                        cur_data = sp_data[bi, :, :]
                        ac.output.nc_write(ofile, b, cur_data, dataset_attributes=att)
                        print('Wrote {} ({})'.format(b, cur_data.shape))
                        cur_data = None
        ## end write Landsat L2 ST

        if not settings['surface_reflectance']:
            ## write thermal bands
            bii = 0
            for bi, b in enumerate(obands_rhot):
                if b not in thermal_bands: continue
                att = {'band': b}
                if rhot_data is not None:
                    cur_data = rhot_data[bi, :, :] * 1.0
                    dso = b.replace('B', 'bt')
                    ac.output.nc_write(ofile, dso, cur_data, dataset_attributes=att)
                    print('Wrote {} ({})'.format(dso, cur_data.shape))
                    cur_data = None

                if st_data is not None:
                    if b in obands_st:
                        if len(st_data.shape) == 3:
                            cur_data = st_data[bii, :, :] * 1.0
                        else:
                            cur_data = st_data * 1.0
                        dso = b.replace('B', 'st')
                        ac.output.nc_write(ofile, dso, cur_data, dataset_attributes=att)
                        print('Wrote {} ({})'.format(dso, cur_data.shape))
                        cur_data = None
                        bii += 1

        print('Wrote {}'.format(ofile))
        rhot_data = None
        rhos_data = None
        sr_data = None
        st_data = None
        geom_data = None

        ## clear intermediate files
        if settings['store_output_locally'] & settings['clear_output_zip_files']:
            for zfile in rhot_files+geom_files+rhos_files+sr_files+sp_files:
                if os.path.exists(zfile):
                    print('Removing {}'.format(zfile))
                    os.remove(zfile)

    return({gatts['acolite_type']:ofile})
