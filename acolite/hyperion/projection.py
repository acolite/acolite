## def get_projection
## get projection from Hyperion metadata
##
## written by Quinten Vanhellemont, RBINS
## 2018-01-31
## modifications: 2021-08-04 (QV) modified to return projection dct and updated for generic acolite

def projection(metadata):
    from pyproj import Proj

    pixel_size = [float(metadata['PROJECTION_PARAMETERS']['GRID_CELL_SIZE'])]
    pixel_size += [pixel_size[0]*-1.0]
    nrow=int(metadata['PRODUCT_METADATA']['PRODUCT_LINES'])
    ncol=int(metadata['PRODUCT_METADATA']['PRODUCT_SAMPLES'])
    dimensions = (ncol,nrow)

    proj = metadata['PROJECTION_PARAMETERS']['MAP_PROJECTION']
    ellipsoid = metadata['PROJECTION_PARAMETERS']['REFERENCE_ELLIPSOID']
    datum = metadata['PROJECTION_PARAMETERS']['REFERENCE_DATUM']

    if (proj == 'UTM') & (ellipsoid == 'WGS84') & (datum == 'WGS84'): is_utm = True
    else: is_utm = False

    if (proj == 'PS') & (ellipsoid == 'WGS84') & (datum == 'WGS84'): is_ps = True
    else: is_ps = False

    if (not is_utm) & (not is_ps): print('Projection not implemented')

    if is_utm:
        zone = int(metadata['UTM_PARAMETERS']['ZONE_NUMBER'])
        south = True if zone < 0 else False
        zone = abs(zone)
        proj4_list = ['+proj=utm',
                      '+zone={}'.format(zone),
                      '+datum={}'.format(datum),
                      '+units=m',
                      '+no_defs ']
        if south: proj4_list += ['+south']

    if is_ps:
        print('Not implemented')

    proj4_string = ' '.join(proj4_list)
    p = Proj(proj4_string)

    ## check corners of image
    x,y = [],[]
    for corner in ['LL','UL','UR','LR']:
            x.append(float(metadata['PRODUCT_METADATA']['PRODUCT_{}_CORNER_MAPX'.format(corner)]))
            y.append(float(metadata['PRODUCT_METADATA']['PRODUCT_{}_CORNER_MAPY'.format(corner)]))
    xrange = [min(x),max(x)]
    yrange = [max(y),min(y)]

    dct = {'p': p, 'epsg':  p.crs.to_epsg(),
           'xrange': xrange, 'yrange': yrange,
           'proj4_string':proj4_string, 'dimensions':dimensions,
           'pixel_size': pixel_size, 'utm': is_utm, 'ps': is_ps}

    if is_utm: dct['zone'] : zone
    if is_ps:
        dct['vertical_lon'] : vertical_lon
        dct['lat_ts'] : lat_ts
        dct['false_e'] : false_e
        dct['false_n'] : false_n
        dct['lat_0'] : lat_0

    dct['xdim'] = int((dct['xrange'][1]-dct['xrange'][0])/dct['pixel_size'][0])+1
    dct['ydim'] = int((dct['yrange'][1]-dct['yrange'][0])/dct['pixel_size'][1])+1

    return(dct)
