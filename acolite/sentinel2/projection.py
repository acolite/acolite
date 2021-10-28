## def get_projection
## get projection from Sentinel granule metadata
##
## written by Quinten Vanhellemont, RBINS
## 2017-04-18
## modifications: 2018-03-07 (QV) changed to Proj4 string to avoid _gdal import problem in PyInstaller
##                2018-06-06 (QV) added return of Proj4 string
##                2021-02-11 (QV) renamed, reformatted for integration in acolite-gen

def projection(meta, s2_target_res=10, return_grids=False):
    from pyproj import Proj

    is_utm = True
    is_ps = False

    proj = 'UTM'
    ellipsoid = 'WGS84'
    datum = 'WGS84'
    cs_code = meta["HORIZONTAL_CS_CODE"]
    cs_name = meta["HORIZONTAL_CS_NAME"]
    epsg = int(cs_code.split(':')[1])

    split = cs_name.split('/')
    datum = split[0].strip()
    zone_name = split[1].split()[-1]

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

    ## construct 10, 20 and 60m grids
    grids = {}
    for i in meta['GRIDS'].keys():
        grid=meta['GRIDS'][i]
        x0 = float(grid['ULX'])
        y0 = float(grid['ULY'])
        xs = float(grid['XDIM'])
        ys = float(grid['YDIM'])
        nx = float(grid['NCOLS'])
        ny = float(grid['NROWS'])

        x1 = x0 + (xs * nx)
        y1 = y0 + (ys * ny)
        xrange = (x0,x1)
        yrange = (y0,y1)
        grids[int(i)] = {'xrange':xrange,'yrange':yrange, 'nx':nx,'ny':ny,
                         'x0':x0, 'xs':xs, 'x1':x1, 'y0':y0, 'ys':ys, 'y1':y1}

    ## select target grid here and return projection dict
    sel_grid = grids[int(s2_target_res)]
    dimensions = sel_grid['nx'], sel_grid['ny']
    pixel_size = sel_grid['xs'], sel_grid['ys']
    dct = {'p': p, 'epsg': p.crs.to_epsg(),
           'xrange': [sel_grid['xrange'][0], sel_grid['xrange'][1]],
           'yrange': [sel_grid['yrange'][0], sel_grid['yrange'][1]],
           'proj4_string':proj4_string, 'dimensions':dimensions,
           'pixel_size': pixel_size, 'utm': is_utm, 'ps': is_ps}
    if is_utm: dct['zone'] : zone_name

    ## offset end of range by one pixel to make correct lat/lon x/y datasets
    #dct['xrange'][1]-=dct['pixel_size'][0]
    #dct['yrange'][1]-=dct['pixel_size'][1]
    #dct['xdim'] = int((dct['xrange'][1]-dct['xrange'][0])/dct['pixel_size'][0])+1
    #dct['ydim'] = int((dct['yrange'][1]-dct['yrange'][0])/dct['pixel_size'][1])+1
    dct['xdim'] = int((dct['xrange'][1]-dct['xrange'][0])/dct['pixel_size'][0])
    dct['ydim'] = int((dct['yrange'][1]-dct['yrange'][0])/dct['pixel_size'][1])

    if return_grids:
        return(dct, grids)
    else:
        return(dct)
