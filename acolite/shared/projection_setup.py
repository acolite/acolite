## def projection_setup
## sets up a suitable utm projection based on SWNE limit
## written by Quinten Vanhellemont, RBINS
## 2022-01-18
## modifications: 2022-07-12 (QV) changed name from projection_utm

def projection_setup(limit, resolution, res_method = 'bilinear', utm=True, epsg=None, add_half_pixel=False):
    from pyproj import Proj
    import acolite as ac
    import numpy as np

    if (utm) and (epsg is None):
        lon = (limit[1]+limit[3])/2
        lat = (limit[0]+limit[2])/2
        utm_zone, epsg = ac.shared.utm_epsg(lon, lat)
    if (~utm) and (epsg is None):
        epsg = 4326

    projection = 'EPSG:{}'.format(epsg)

    p = Proj(projection)
    if len(np.atleast_1d(resolution)) == 1:
        target_pixel_size = resolution, resolution * -1.0
    else:
        target_pixel_size = resolution

    xrange_raw, yrange_raw = p((limit[1],limit[1],limit[3],limit[3]),
                               (limit[0],limit[2],limit[2],limit[0]))
    xrange_raw = (min(xrange_raw), max(xrange_raw))
    yrange_raw = (min(yrange_raw), max(yrange_raw))

    xrange_region = [xrange_raw[0] - (xrange_raw[0] % target_pixel_size[0]), \
                     xrange_raw[1]+target_pixel_size[0]-(xrange_raw[1] % target_pixel_size[0])]
    yrange_region = [yrange_raw[1]+target_pixel_size[1]-(yrange_raw[1] % target_pixel_size[1]), \
                     yrange_raw[0] - (yrange_raw[0] % target_pixel_size[1])]

    x_grid_off = xrange_region[0]%target_pixel_size[0], xrange_region[1]%target_pixel_size[0]
    y_grid_off = yrange_region[0]%target_pixel_size[1], yrange_region[1]%target_pixel_size[1]
    xrange = (xrange_region[0]-x_grid_off[0]), (xrange_region[1]+(target_pixel_size[0]-x_grid_off[1]))
    yrange = (yrange_region[0]-y_grid_off[0]), (yrange_region[1]+(target_pixel_size[1]-y_grid_off[1]))

    if add_half_pixel:
        xrange = xrange[0]-(target_pixel_size[0]/2), xrange[1]-(target_pixel_size[0]/2)
        yrange = yrange[0]-(target_pixel_size[1]/2), yrange[1]-(target_pixel_size[1]/2)

    ## pixel sizes
    #ny = int((yrange[0] - yrange[1])/target_pixel_size[1])
    nx = int(np.round((xrange[1] - xrange[0])/target_pixel_size[0]))
    ny = int(np.round((yrange[1] - yrange[0])/target_pixel_size[1]))

    ## set up projection dict and nc_projection
    dct = {'xrange': xrange, 'yrange': yrange, 'p': p,
           'pixel_size': target_pixel_size, 'xdim': nx, 'ydim': ny, 'epsg':projection, 'projection':projection}
    nc_projection = ac.shared.projection_netcdf(dct, add_half_pixel=True)

    xyr = [min(dct['xrange']),
           min(dct['yrange']),
           max(dct['xrange']),
           max(dct['yrange']),
           projection]

    warp_to = (projection, xyr, dct['pixel_size'][0], dct['pixel_size'][1], res_method)
    return(dct, nc_projection, warp_to)
