## def projection_utm
## sets up a suitable utm projection based on SWNE limit
## written by Quinten Vanhellemont, RBINS
## 2022-01-18
## modifications:

def projection_utm(limit, resolution, res_method = 'bilinear'):
    from pyproj import Proj
    import acolite as ac

    lon = (limit[1]+limit[3])/2
    lat = (limit[0]+limit[2])/2
    utm_zone, epsg = ac.shared.utm_epsg(lon, lat)

    projection = 'EPSG:{}'.format(epsg)

    p = Proj(projection)
    target_pixel_size = resolution, resolution * -1.0

    xrange_raw, yrange_raw = p((limit[1],limit[1],limit[3],limit[3]),
                               (limit[0],limit[2],limit[2],limit[0]))
    xrange_raw = (min(xrange_raw), max(xrange_raw))
    yrange_raw = (min(yrange_raw), max(yrange_raw))
    xrange_region = [xrange_raw[0] - (xrange_raw[0] % target_pixel_size[0]*2), xrange_raw[1]+target_pixel_size[0]*2-(xrange_raw[1] % target_pixel_size[0]*2)]
    yrange_region = [yrange_raw[1]+target_pixel_size[1]*2-(yrange_raw[1] % target_pixel_size[1]*2), yrange_raw[0] - (yrange_raw[0] % target_pixel_size[1]*2)]


    x_grid_off = xrange_region[0]%target_pixel_size[0], xrange_region[1]%target_pixel_size[0]
    y_grid_off = yrange_region[0]%target_pixel_size[1], yrange_region[1]%target_pixel_size[1]
    xrange = (xrange_region[0]-x_grid_off[0]), (xrange_region[1]+(target_pixel_size[0]-x_grid_off[1]))
    yrange = (yrange_region[0]-y_grid_off[0]), (yrange_region[1]+(target_pixel_size[0]-y_grid_off[1]))

    ## pixel sizes
    ny = int((yrange[0] - yrange[1])/target_pixel_size[1])
    nx = int((xrange[1] - xrange[0])/target_pixel_size[0])
    ny = int((yrange[1] - yrange[0])/target_pixel_size[1])

    ## set up projection dict and nc_projection
    dct = {'xrange': xrange, 'yrange': yrange, 'p': p,
           'pixel_size': target_pixel_size, 'xdim': nx, 'ydim': ny}
    nc_projection = ac.shared.projection_netcdf(dct, add_half_pixel=True)

    xyr = [min(dct['xrange']),
           min(dct['yrange']),
           max(dct['xrange']),
           max(dct['yrange']),
           projection]

    warp_to = (projection, xyr, dct['pixel_size'][0], dct['pixel_size'][1], res_method)
    return(dct, nc_projection, warp_to)
