## def l1_convert
## reads radiance metadata and strip_id from Tanager L1B RAD file
## written by Quinten Vanhellemont, RBINS
## 2025-07-22
## modifications: 

def metadata(bundle):
    from netCDF4 import Dataset
    meta = {}
    with Dataset(bundle) as nc:
        file_attrs = nc.groups['HDFEOS']['SWATHS']['HYP'].ncattrs()
        for a in file_attrs:
            meta[a] = nc.groups['HDFEOS']['SWATHS']['HYP'].getncattr(a)

        ## get radiance attributes
        rad_attrs = nc.groups['HDFEOS']['SWATHS']['HYP']['Data Fields']['toa_radiance'].ncattrs()
        for a in rad_attrs:
            meta[a] = nc.groups['HDFEOS']['SWATHS']['HYP']['Data Fields']['toa_radiance'].getncattr(a)
    return(meta)
