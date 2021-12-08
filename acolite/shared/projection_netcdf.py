## def projection_netcdf
## generates projection dict to use in NetCDF writing
## written by Quinten Vanhellemont, RBINS
## 2021-12-08
## modifications:

def projection_netcdf(dct, add_half_pixel = False):
    import numpy as np

    cf = dct['p'].crs.to_cf()
    proj_key = cf['grid_mapping_name']

    ## nc projection dict
    nc_projection = {'x': {}, 'y': {}, proj_key: {'data':None, 'attributes': cf}}

    ## add x and y
    nc_projection['x']['attributes'] = {'standard_name': 'projection_x_coordinate',
                                        'long_name': 'x coordinate of projection', 'units': 'm'}
    nc_projection['y']['attributes'] = {'standard_name': 'projection_y_coordinate',
                                        'long_name': 'y coordinate of projection', 'units': 'm'}

    ## compute x and y vectors
    if not add_half_pixel:
        x = np.linspace(dct['xrange'][0], dct['xrange'][1]-dct['pixel_size'][0], dct['xdim'])
        y = np.linspace(dct['yrange'][0], dct['yrange'][1]-dct['pixel_size'][1], dct['ydim'])
    else:
        x = np.linspace(dct['xrange'][0],dct['xrange'][1]-dct['pixel_size'][0],dct['xdim']) + dct['pixel_size'][0]/2
        y = np.linspace(dct['yrange'][0],dct['yrange'][1]-dct['pixel_size'][1],dct['ydim']) + dct['pixel_size'][1]/2
    nc_projection['x']['data'] = x
    nc_projection['y']['data'] = y

    return(nc_projection)
