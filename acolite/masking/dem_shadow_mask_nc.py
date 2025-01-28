## def dem_shadow_mask_nc
## makes shadow mask for an ACOLITE NetCDF file
## based on the Richens 1997 shadow volume algorithm "Image processing for urban scale environmental modelling, Paul Richens, 1997b"
## written by Quinten Vanhellemont, RBINS
## 2024-01-25
## modifications: 2024-01-26 (QV) added option to extend grid in the direction of the sun
##                                updated determination of pixel_size, added acolite settings
##                2024-01-31 (QV) updated grid extension, allow for 4 element extension to be provided [top, bottom, left, right]
##                2024-02-01 (QV) check DEM size, check binary operations iteration

def dem_shadow_mask_nc(ncf, return_dem = False, extend = []):
    import acolite as ac
    import numpy as np
    import dateutil.parser, pytz, scipy.ndimage

    if ac.settings['run']['verbosity'] > 2: print('Running DEM cast shadow mask for {}'.format(ncf))

    ## read datasets and gatts
    datasets = ac.shared.nc_datasets(ncf)
    gatts = ac.shared.nc_gatts(ncf)

    ## read lat/lon
    lat = ac.shared.nc_data(ncf, 'lat')
    lon = ac.shared.nc_data(ncf, 'lon')

    ## compute scene centre lon/lat
    centre_lon = np.nanmedian(lon)
    centre_lat = np.nanmedian(lat)

    ## get nc projection
    nc_projection = ac.shared.nc_projection_read(ncf)
    if (nc_projection is None) & (ac.settings['run']['dem_shadow_mask_extend']):
        if ac.settings['run']['verbosity'] > 2: print('Could not read nc_projection, not extending grid.')
        ac.settings['run']['dem_shadow_mask_extend'] = False

    ## get pixel size
    if nc_projection is not None: ## from nc projection
        dct = ac.shared.nc_projection_dct(nc_projection)
        pixel_size = 1.0 * dct['pixel_size'][0]
        if ac.settings['run']['verbosity'] > 2: print('Using pixel_size={} for sensor {} from nc_projection'.format(pixel_size, gatts['sensor']))
    elif 'pixel_size' in gatts: ## from gatts
        pixel_size = gatts['pixel_size'][0]
        if ac.settings['run']['verbosity'] > 2: print('Using pixel_size={} for sensor {} from attributes'.format(pixel_size, gatts['sensor']))
    else: ## assume default pixel sizes
        if ('OLI' in gatts['sensor']) | ('TM' in gatts['sensor']) | ('ETM' in gatts['sensor']):
            pixel_size = 30
        elif 'MSI' in gatts['sensor']:
            pixel_size = 10
        elif 'PlanetScope' in gatts['sensor']:
            pixel_size = 3
        else:
            if ac.settings['run']['verbosity'] > 0: print('Could not determine pixel_size for sensor {}'.format(gatts['sensor']))
            return

    ## get scene date
    isodate = gatts['isodate']
    date = dateutil.parser.parse(isodate)
    date = date.replace(tzinfo=pytz.UTC)

    ## get sun position
    spos = ac.shared.sun_position(date, centre_lon, centre_lat)
    sza = spos['zenith'][0]
    saa = spos['azimuth'][0]

    ## estimate grid convergence to adjust sun azimuth
    zone = None
    gca = 0
    if 'mgrs_tile' in gatts: ## Sentinel-2 tiling grid SAFE
        zone = int(gatts['mgrs_tile'][1:3])
    elif 'MGRS_TILE' in gatts: ## Sentinel-2 tiling grid GEE
        zone = int(gatts['MGRS_TILE'][0:2])
    elif 'scene_zone' in gatts: ## Landsat zone
        zone = int(gatts['scene_zone'])
    elif 'proj4_string' in gatts:
        if 'stere' in gatts['proj4_string']:
            zone = None
        elif 'utm' in gatts['proj4_string']:
            zone = int([k.split('=')[1] for k in gatts['proj4_string'].split() if 'zone' in k][0])
    else:
        ## guess zone?
        zone = int(ac.shared.utm_epsg(centre_lon, centre_lat)[0])
        if ac.settings['run']['verbosity'] > 2: print('Estimated zone {} from scene centre {:.3f}°N {:.3f}°E'.format(zone, centre_lat, centre_lon))

    ## get grid convergence
    if zone is not None:
        ## central meridian
        lon0 = (6 * zone - 183)
    else:
        lon0 = 0 ## at least for southern hemisphere stere?

    # ## estimate projection azimuth
    # lat1, lon1 = lat[-1, int(lat.shape[1]/2)], lon[-1, int(lat.shape[1]/2)]
    # lat2, lon2 = lat[0, int(lat.shape[1]/2)], lon[0, int(lat.shape[1]/2)]
    # paa = ac.shared.azimuth_two_points(lon1, lat1, lon2, lat2)
    # if paa > 180: paa = paa-360
    # print(paa)

    ## grid convergence
    gc = np.arctan(np.tan(np.radians(centre_lon) - np.radians(lon0)) * np.sin(np.radians(centre_lat)))
    gca = np.degrees(gc)

    ## correct saa for grid convergence
    saa_ = saa-gca
    if saa_ < 0: saa_ += 360

    if ac.settings['run']['verbosity'] > 2:
        print('Grid convergence = {:.2f} degrees'.format(gca))
        print('Sun azimuth = {:.2f} degrees'.format(saa))
        print('Sun azimuth adjusted for grid convergence = {:.2f} degrees'.format(saa_))

    ## get DEM
    if ('dem' in datasets) & (ac.settings['run']['dem_shadow_mask_extend'] is False):
        if ac.settings['run']['verbosity'] > 2: print('Reading DEM data from {}'.format(ncf))
        dem = ac.shared.nc_data(ncf, 'dem')
    else:
        if ac.settings['run']['verbosity'] > 2: print('Getting DEM data from {}'.format(ac.settings['run']['dem_source']))
        if ac.settings['run']['dem_shadow_mask_extend']:
            ## amount of pixels to extend the image
            xext = int(np.round(ac.settings['run']['dem_shadow_mask_extend_metres']/dct['pixel_size'][0]))
            yext = int(np.round(ac.settings['run']['dem_shadow_mask_extend_metres']/dct['pixel_size'][1]))

            ## find directions to extend
            extend_top = False
            extend_bottom = False
            extend_left = False
            extend_right = False
            if (saa_ >= 0) & (saa_ < 90):
                extend_top = True
                extend_right = True
            if (saa_ >= 90) & (saa_ < 180):
                extend_bottom = True
                extend_right = True
            if (saa_ >= 180) & (saa_ < 270):
                extend_bottom = True
                extend_left = True
            if (saa_ >= 270) & (saa_ < 360):
                extend_top = True
                extend_left = True

            ## if requested extension
            if len(extend) == 4:
                extend_top, extend_bottom, extend_left, extend_right = extend

            ## get extended ranges
            xrange_ = [dct['xrange'][0], dct['xrange'][1]]
            yrange_ = [dct['yrange'][0], dct['yrange'][1]]

            if extend_right:
                xrange_[1] = dct['xrange'][1] + dct['pixel_size'][0] * xext
            if extend_left:
                xrange_[0] = dct['xrange'][0] - dct['pixel_size'][0] * xext
            if extend_top:
                yrange_[0] = dct['yrange'][0] + dct['pixel_size'][1] * yext
            if extend_bottom:
                yrange_[1] = dct['yrange'][1] - dct['pixel_size'][1] * yext

            ## extended dimensions
            xdim_ = int((xrange_[1]-xrange_[0]+dct['pixel_size'][0])/dct['pixel_size'][0])
            ydim_ = int((yrange_[1]-yrange_[0]+dct['pixel_size'][1])/dct['pixel_size'][1])

            ## subsetting for extended array
            y0 = int(np.abs((dct['yrange'][0]-yrange_[0])/dct['pixel_size'][1]))
            y1 = int(ydim_ - np.abs((dct['yrange'][1]-yrange_[1])/dct['pixel_size'][1]))
            x0 = int(np.abs((dct['xrange'][0]-xrange_[0])/dct['pixel_size'][0]))
            x1 = int(xdim_ - np.abs((dct['xrange'][1]-xrange_[1])/dct['pixel_size'][0]))

            if ac.settings['run']['verbosity'] > 3:
                print('Extending top={} bottom={} left={} right={}'\
                        .format(extend_top, extend_bottom, extend_left, extend_right))
                print('yext={} xext={}'.format(yext, xext))
                print('ydim_={} xdim_={}'.format(ydim_, xdim_))

            ## update dct
            dct_ = {k:dct[k] for k in dct}
            dct_['xrange'] = xrange_
            dct_['yrange'] = yrange_
            dct_['xdim'] = xdim_
            dct_['ydim'] = ydim_
            dct_['dimensions'] = (xdim_, ydim_)

            ## get extended geolocation
            lon_, lat_ = ac.shared.projection_geo(dct_)
            dem = ac.dem.dem_lonlat(lon_, lat_, source=ac.settings['run']['dem_source'])

            ## test dem shape (e.g. SRTM is limited in latitude)
            if dem.shape == ():
                print('Could not retrieve DEM from source {}'.format(ac.settings['run']['dem_source']))
                return
        else:
            dem = ac.dem.dem_lonlat(lon, lat, source=ac.settings['run']['dem_source'])

    ## compute shadow mask
    if ac.settings['run']['verbosity'] > 2: print('Computing terrain cast shadows')
    mask_shadow = ac.masking.dem_shadow_mask(dem, saa_, sza, pixel_size)
    shade = np.ma.masked_where(mask_shadow == 0, mask_shadow)

    ## dilate erode mask?
    if ac.settings['run']['dem_shadow_mask_filter']:
        if ac.settings['run']['dem_shadow_mask_erode_its'] > 0:
            mask_shadow = scipy.ndimage.binary_erosion(mask_shadow, np.ones(ac.settings['run']['dem_shadow_mask_filter_kernel']),
                                                                     iterations = ac.settings['run']['dem_shadow_mask_erode_its'])

        if ac.settings['run']['dem_shadow_mask_dilate_its'] > 0:
            mask_shadow = scipy.ndimage.binary_dilation(mask_shadow, np.ones(ac.settings['run']['dem_shadow_mask_filter_kernel']),
                                                                 iterations = ac.settings['run']['dem_shadow_mask_dilate_its'])

        ## smooth mask?
        #mask_shadow = scipy.ndimage.gaussian_filter(mask_shadow, 1, mode='reflect')

        ## mask mask
        shade = np.ma.masked_where(mask_shadow == 0, mask_shadow)

    ## return dem
    if (return_dem):
        if ac.settings['run']['dem_shadow_mask_extend']:
            return(shade, dem, lon_, lat_, y0, y1, x0, x1)
        else:
            return(shade, dem, lon, lat)

    ## crop if image was extended
    if ac.settings['run']['dem_shadow_mask_extend']: shade = shade[y0:y1,x0:x1]

    return(shade)
