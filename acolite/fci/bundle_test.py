## def bundle_test
## test if given bundle has FCI file(s)
## basically read the first nc file and check if the data source is FCI
##
## written by Quinten Vanhellemont, RBINS
## 2025-05-13
## modifications:

def bundle_test(bundle_in, match = '*.nc'):
    import os, glob
    import acolite as ac

    if os.path.isdir(bundle_in):
        bundle_dir = bundle_in
        fn = None
    else:
        bundle_dir = os.path.dirname(bundle_in)
        fn = os.path.basename(bundle_in)

    ## list files
    files = glob.glob('{}/{}'.format(bundle_dir, match))
    files.sort()

    sel_nc = None
    for file in files:
        if (os.path.splitext(file)[-1] == '.nc'):
            if fn is None:
                fn = os.path.basename(file)
                sel_nc = '{}'.format(file)
            else:
                if fn !=  os.path.basename(file): continue
                sel_nc = '{}'.format(file)
            break

    if sel_nc is not None:
        gatts = ac.shared.nc_gatts(sel_nc)
        if gatts['data_source'] == 'FCI':
            return(sel_nc)

    return
