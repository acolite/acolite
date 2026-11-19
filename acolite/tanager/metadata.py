## def l1_convert
## reads radiance metadata and strip_id from Tanager L1B RAD file
## written by Quinten Vanhellemont, RBINS
## 2025-07-22
## modifications: 2026-11-19 (QV) added support for ortho_radiance

def metadata(bundle):
    from netCDF4 import Dataset
    meta = {}
    with Dataset(bundle) as nc:
        if 'SWATHS' in nc['HDFEOS'].groups.keys(): group_key = 'SWATHS'
        if 'GRIDS' in nc['HDFEOS'].groups.keys(): group_key = 'GRIDS'

        file_attrs = nc.groups['HDFEOS'][group_key]['HYP'].ncattrs()
        for a in file_attrs:
            meta[a] = nc.groups['HDFEOS'][group_key]['HYP'].getncattr(a)

        ## get radiance attributes
        rad_attrs = nc.groups['HDFEOS'][group_key]['HYP']['Data Fields']['toa_radiance'].ncattrs()
        for a in rad_attrs:
            meta[a] = nc.groups['HDFEOS'][group_key]['HYP']['Data Fields']['toa_radiance'].getncattr(a)
    return(meta)
