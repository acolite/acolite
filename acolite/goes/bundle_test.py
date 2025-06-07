## def bundle_test
## test if given bundle has ABI file(s)
##
## written by Quinten Vanhellemont, RBINS
## 2025-06-05
## modifications:

def bundle_test(bundle, rad_base = 'OR_ABI-L1b-Rad'):
    import os, glob
    import acolite as ac

    if os.path.isdir(bundle):
        files = glob.glob('{}/{}*-*.nc'.format(bundle, rad_base))
    else:
        bn = os.path.basename(bundle)
        dn = os.path.dirname(bundle)
        sp = bn.split('_')
        files = glob.glob('{}/{}*-*_{}_{}_*.nc'.format(dn, rad_base, sp[2], sp[3]))
    files.sort()

    bundle_files = {}

    for file in files:
        bn = os.path.basename(file)
        sp = bn.split('_')
        satellite = sp[2]
        start = sp[3]
        sp = sp[1].split('-')
        target = sp[2]
        band = sp[3]
        #print(satellite, target, start, band)
        if satellite not in bundle_files: bundle_files[satellite] = {}
        if target not in bundle_files[satellite]: bundle_files[satellite][target] = {}
        if start not in bundle_files[satellite][target]: bundle_files[satellite][target][start] = {}
        if band in bundle_files[satellite][target][start]:
            print('Band {} already in bundle files'.format(band))
            print('Old path {}'.format(bundle_files[satellite][target][start][band]['path']))
            print('New path {}'.format(file))

        bundle_files[satellite][target][start][band] = {'path': file, 'gatts': ac.shared.nc_gatts(file)}

    return(bundle_files)
