## def projection
## get projection dct from Tanager ortho_radiance bundle
## written by Quinten Vanhellemont, RBINS
## 2026-11-19
## modifications:

def projection(bundle):
    import netCDF4
    from pyproj import Proj

    with netCDF4.Dataset(bundle, 'r') as nc:
        struct_string = nc['HDFEOS INFORMATION']['StructMetadata.0'][()]
        epsg_code = nc['HDFEOS']['GRIDS']['HYP'].getncattr('epsg_code')

    projection = 'EPSG:{}'.format(epsg_code)
    p = Proj(projection)


    ## parse the metadata
    struct_list = []
    for line in struct_string.split('\n'):
        line = line.strip('\t')
        struct_list.append(line)
    struct_dict = {}
    for l in struct_list:
        sp = l.split('=')
        if len(sp) != 2: continue
        if sp[0] in ['GROUP', 'END_GROUP', 'OBJECT', 'END_OBJECT']:
            continue
        if sp[0] not in struct_dict:
            if sp[0] in ['UpperLeftPointMtrs', 'LowerRightMtrs', 'DimList', 'MaxdimList']:
                struct_dict[sp[0]] = sp[1][1:-2].split(',')
            else:
                struct_dict[sp[0]] = sp[1]

    xrange = [float(struct_dict['UpperLeftPointMtrs'][0]), float(struct_dict['LowerRightMtrs'][0])]
    yrange = [float(struct_dict['UpperLeftPointMtrs'][1]), float(struct_dict['LowerRightMtrs'][1])]

    nx = int(struct_dict['XDim'])
    ny = int(struct_dict['YDim'])
    target_pixel_size = (xrange[1]-xrange[0])/nx, (yrange[1]-yrange[0])/ny

    ## set up projection dict and nc_projection
    dct = {'xrange': xrange, 'yrange': yrange, 'p': p,
           'pixel_size': target_pixel_size, 'xdim': nx, 'ydim': ny, 'epsg':projection, 'projection':projection}

    return(dct)
