## Helper for converting WISE data to ACOLITE l1R NetCDF
## written by Raphael Mabit, ISMER-UQAR
## 2023-12-12

def rm_illegal_chr_envi_hdr(header):
    '''
    remove illegal characters in Envi header, which cause spectral's FileNotEnviHeader exception
    :param header:
    :return:
    '''
    with open(header, 'r', encoding='latin-1') as f:
        lines = f.readlines()
        for i, l in enumerate(lines):
            if l.find("²/")  > -1:
                lines[i] = l.replace("²/", '')
            if l.find("± ") > -1:
                lines[i] = l.replace("± ", '')
        return lines

def parse_mapinfo(header):
    '''
    compose affine from map info in ENVI Header
    :param:  map info, a list like [UTM, 1, 1, 552591.000, 5450146.500, 1.500, 1.500, 19, North,  WGS84]
    :param:
    :return:  Affine  (from rasterio.transform import Affine)
    '''
    import pyproj

    proj_name, x, y, reference_x, reference_y, x_size, y_size, proj_zone, north_or_south, datum = header['map info'].replace(" ", "").split(',')

    ## get projection info
    # PROJ string need the flag +south if utm in the southern hemisphere
    # https://proj.org/en/9.3/operations/projections/utm.html
    north_or_south = north_or_south[0]
    unit = 'm'
    proj4string  = "+proj={0} +zone={1} +ellps={2} +datum={2} +units={3} +no_defs".format(
        str.lower(proj_name),
        proj_zone,
        str.upper(datum),
        unit)

    crs = pyproj.CRS.from_proj4(proj4string)
    # The CRS doesn't have all the same proprieties if created from PROJ string or EPSG
    # You can notice that in the conversion to WTK
    crs = pyproj.CRS.from_epsg(crs.to_epsg())
    Wkt = crs.to_wkt()
    p = pyproj.Proj(Wkt)

    #### FOR DEV
    # Code taken from shared.projection_read.py
    # TODO check pixel reference, upper left or center ?
    #   X and Y coordinates of origin pixel are given as
    #   upper left corner of upper left pixel by GDAL (https://gdal.org/user/raster_data_model.html)
    #   Seems to have an incoherence between the pixel center of wtransform and trasform 4 5
    # if wtransform is not None:
    #     ## world transform elements:
    #     ## 0 pixel X size
    #     ## 1 rotation about the Y axis (usually 0.0)
    #     ## 2 rotation about the X axis (usually 0.0)
    #     ## 3 negative pixel Y size
    #     ## 4 X coordinate of upper left pixel center
    #     ## 5 Y coordinate of upper left pixel center
    #     x0 = wtransform[4]
    #     dx = wtransform[0]
    #     y0 = wtransform[5]
    #     dy = wtransform[3]
    # else:
    #     ## derive projection extent

    #     x0 = transform[0]
    #     dx = transform[1]
    #     y0 = transform[3]
    #     dy = transform[5]
    #### END DEV

    # Raster Size
    dimx = int(header['samples']) # Rows
    dimy = int(header['lines']) # Columns

    # Pixel origin, upper left pixel
    # ENVI header map info description doesn't tell about the tie point location (Upper Left ? Center ?)
    # https://www.nv5geospatialsoftware.com/docs/ENVIHeaderFiles.html
    # In ACWISE I found the comment:
    # "1 1 means the referring point is (left,upper), this is the most often case"
    x0 = float(reference_x)
    y0 = float(reference_y)

    # Pixel size
    dx = float(x_size)
    dy = -float(y_size)

    pixel_size = (dx, dy)

    xrange = (x0, x0+dimx*dx)
    yrange = (y0, y0+dimy*dy)

    ## make acolite generic dict
    dct = {'p': p, 'epsg': p.crs.to_epsg(),
           'Wkt': Wkt,  'proj4_string': crs.to_proj4(),
           'xrange': xrange, 'yrange': yrange,
           'xdim':dimx, 'ydim': dimy,
           'dimensions':(dimx, dimy),
           'pixel_size': pixel_size}
    dct['projection'] = 'EPSG:{}'.format(dct['epsg'])

    return(dct)
