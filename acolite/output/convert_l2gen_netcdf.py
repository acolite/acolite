## convert_l2gen_netcdf
## function to convert l2gen NetCDF to something useable with ACOLITE functions
## (i.e. remove geophysical and navigation data groups)
##
## written by Quinten Vanhellemont, RBINS
## 2023-05-31
## modifications:

def convert_l2gen_netcdf(ifile, output=None, file_type = 'L2S', \
                         l2_mask_value = 0, rotate = False, reproject = False, settings = {}):
    import os
    import numpy as np
    import acolite as ac

    ## output paths
    dn = os.path.dirname(ifile)
    bn = os.path.basename(ifile)
    bno = bn.replace('.nc', '_{}.nc'.format(file_type))
    ofile = '{}/{}'.format(output if output is not None else dn, bno)
    if not os.path.exists(os.path.dirname(ofile)): os.makedirs(os.path.dirname(ofile))

    ## read gatts and dataset info
    gatts = ac.shared.nc_gatts(ifile)
    if gatts['publisher_name'] != 'NASA/GSFC/OBPG':
        print('Publisher of {} not NASA/GSFC/OBPG')
        print('Exiting.')
        return
    datasets = ac.shared.nc_datasets(ifile, group='geophysical_data')
    datasets_nav = ac.shared.nc_datasets(ifile, group='navigation_data')

    ## add sensor for ACOLITE funs
    gatts['sensor'] = '{}_{}'.format(gatts['platform'], gatts['instrument']).upper()

    ## determine additional masking
    l2_flags, att = ac.shared.nc_data(ifile, 'l2_flags', group='geophysical_data', attributes=True)
    l2_mask = (l2_flags & l2_mask_value) != 0

    new = True
    for ds in datasets_nav:
        if ds not in ['longitude', 'latitude']: continue
        d, att = ac.shared.nc_data(ifile, ds, group='navigation_data', attributes=True)
        if rotate: d = np.rot90(d, k=2)
        ac.output.nc_write(ofile, ds[0:3], d, dataset_attributes=att, new=new, attributes=gatts)
        new = False
        print('Wrote {} to {}'.format(ds[0:3], ofile))

    for ds in datasets:
        if ds in ['l2_flags']: continue
        d, att = ac.shared.nc_data(ifile, ds, group='geophysical_data', attributes=True)
        ## apply mask
        if d.dtype == 'float32':
            d[d.mask] = np.nan
            d[l2_mask] = np.nan
            d = d.data
        else:
            d = d.data

        if rotate: d = np.rot90(d, k=2)
        ac.output.nc_write(ofile, ds, d, dataset_attributes=att, new=new, attributes=gatts)
        new = False
        print('Wrote {} to {}'.format(ds, ofile))

    ## reproject data if requested
    if reproject:
        ac.output.project_acolite_netcdf(ofile, output = None, settings = settings)

    return(ofile)
