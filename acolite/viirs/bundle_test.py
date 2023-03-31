## def bundle_test
## identifies files for VIIRS processing
## written by Quinten Vanhellemont, RBINS
## 2023-03-30
## modifications:

def bundle_test(bundle):
    import os, glob
    import acolite as ac

    bn = os.path.basename(bundle)
    dn = os.path.dirname(bundle)
    meta = ac.shared.nc_gatts(bundle)

    files = None
    if 'VIIRS' in meta['title']:
        sp = bn.split('.')

        ss = sp[0][0:3]
        sl = sp[0][3:5]
        si = sp[0][5:]

        files = {}
        if ss in ['VNP', 'VJ1', 'VJ2']:
            ftypes = {'l1b_mod': '02MOD','geo_mod': '03MOD',
                      'l1b_img': '02IMG','geo_img': '03IMG'}

            for ft in ftypes:
                fs = glob.glob('{}/{}*.nc'.format(dn, '{}*{}.'.format(ss, ftypes[ft])+'.'.join(bn.split('.')[1:4])))
                if len(fs) == 1: files[ft] = fs[0]
    return(files)
