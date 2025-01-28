## def nc_projection_read
## reads projection from netcdf if present
## written by Quinten Vanhellemont, RBINS
## 2021-12-08
## modifications: 2025-01-25 (QV) renamed from nc_read_projection

def nc_projection_read(ncf):
    import acolite as ac
    gatts = ac.shared.nc_gatts(ncf)
    if 'projection_key' not in gatts:
        return(None)
    pkey = gatts['projection_key']
    p, patt = ac.shared.nc_data(ncf, pkey, attributes=True)
    x, xatt = ac.shared.nc_data(ncf, 'x', attributes=True)
    y, yatt = ac.shared.nc_data(ncf, 'y', attributes=True)
    nc_projection = {'x' : {'data': x, 'attributes': xatt},
                     'y' : {'data': y, 'attributes': yatt},
                     pkey : {'data': p, 'attributes': patt}}
    return(nc_projection)
