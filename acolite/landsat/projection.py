## def projection
## gets projection from landsat metadata
## written by Quinten Vanhellemont, RBINS
## 2021-02-05
## modifications: 2021-02-05 (QV) now xrange/yrange represents UL->LR, with negative y pixel size
##                2021-02-05 (QV) added to_epsg

def projection(meta):
    from pyproj import Proj

    if 'PRODUCT_CONTENTS' in meta: # COLL2
        prk = 'PROJECTION_ATTRIBUTES'
        pk = 'PROJECTION_ATTRIBUTES'
    elif 'PRODUCT_METADATA' in meta: ## COLL1
        prk = 'PROJECTION_PARAMETERS'
        pk = 'PRODUCT_METADATA'

    dimensions = int(meta[pk]['REFLECTIVE_LINES']), int(meta[pk]['REFLECTIVE_SAMPLES']) ## Y, X
    pixel = meta[prk]['GRID_CELL_SIZE_REFLECTIVE'] if 'GRID_CELL_SIZE_REFLECTIVE' in meta[prk] \
            else meta[prk]['GRID_CELL_SIZE_REFL']
    pixel = float(pixel)
    pixel_size = pixel, -1 * pixel
    proj = meta[prk]['MAP_PROJECTION']

    ellipsoid = meta[prk]['ELLIPSOID'] if 'ELLIPSOID' in meta[prk] else meta[prk]['REFERENCE_ELLIPSOID']
    datum = meta[prk]['DATUM'] if 'DATUM' in meta[prk] else meta[prk]['REFERENCE_DATUM']

    if (proj == 'UTM') & (ellipsoid == 'WGS84') & (datum == 'WGS84'): is_utm = True
    else: is_utm = False

    if (proj == 'PS') & (ellipsoid == 'WGS84') & (datum == 'WGS84'): is_ps = True
    else: is_ps = False

    if (not is_utm) & (not is_ps): print('Projection not implemented')

    if is_utm:
        zone = int(meta[prk]['UTM_ZONE'])

        #p = Proj(proj='utm',zone=zone,ellps=ellipsoid)
        proj4_list = ['+proj=utm',
                      '+zone={}'.format(abs(zone)),
                      '+datum={}'.format(datum),
                      '+units=m',
                      '+no_defs ']
        if zone < 0: proj4_list += ['+south']


    if is_ps:
        vertical_lon = float(meta[prk]["VERTICAL_LON_FROM_POLE"])
        lat_ts = float(meta[prk]["TRUE_SCALE_LAT"])
        false_e = int(meta[prk]["FALSE_EASTING"])
        false_n = int(meta[prk]["FALSE_NORTHING"])

        lat_0 = -90. if lat_ts < 0 else 90.
        proj4_list =['+proj=stere',
                     '+lat_0={}'.format(lat_0),
                     '+lat_ts={}'.format(lat_ts),
                     '+lon_0={}'.format(vertical_lon),
                     '+k=1',
                     '+x_0={}'.format(false_e),
                     '+y_0={}'.format(false_n),
                     '+datum={}'.format(datum),
                     '+units=m',
                     '+no_defs ']

    proj4_string = ' '.join(proj4_list)
    p = Proj(proj4_string)

    ## check corners of image
    x,y = [],[]
    for corner in ['LL','UL','UR','LR']:
            xtag = 'CORNER_{}_PROJECTION_X_PRODUCT'.format(corner)
            if xtag in meta[pk]:
                x.append(float(meta[pk][xtag]))
            else:
                xtag = 'PRODUCT_{}_CORNER_MAPX'.format(corner)
                x.append(float(meta[pk][xtag]))
            ytag = 'CORNER_{}_PROJECTION_Y_PRODUCT'.format(corner)
            if ytag in meta[pk]:
                y.append(float(meta[pk][ytag]))
            else:
                ytag = 'PRODUCT_{}_CORNER_MAPY'.format(corner)
                y.append(float(meta[pk][ytag]))

    xrange = [min(x)-pixel_size[0]/2,max(x)+pixel_size[0]/2]
    yrange = [max(y)-pixel_size[1]/2,min(y)+pixel_size[1]/2]


    dct = {'p': p, 'epsg':  p.crs.to_epsg(),
           'xrange': xrange, 'yrange': yrange,
           'proj4_string':proj4_string, 'dimensions':dimensions,
           'pixel_size': pixel_size, 'utm': is_utm, 'ps': is_ps}
    if is_utm:
        dct['zone'] = zone
    elif is_ps:
        dct['vertical_lon'] = vertical_lon
        dct['lat_ts'] = lat_ts
        dct['false_e'] = false_e
        dct['false_n'] = false_n
        dct['lat_0'] = lat_0

    dct['xdim'] = int((dct['xrange'][1]-dct['xrange'][0])/dct['pixel_size'][0])#+1
    dct['ydim'] = int((dct['yrange'][1]-dct['yrange'][0])/dct['pixel_size'][1])#+1

    return(dct)
