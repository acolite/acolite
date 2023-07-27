## def bundle_test
## identifies files for VIIRS processing
## written by Quinten Vanhellemont, RBINS
## 2023-03-30
## modifications: 2023-03-31 (QV) support for ocssw processed data/filenames

def bundle_test(bundle):
    import os, glob
    import acolite as ac

    bn = os.path.basename(bundle)
    dn = os.path.dirname(bundle)
    meta = ac.shared.nc_gatts(bundle)

    files = None
    if 'VIIRS' in meta['title']:
        sp = bn.split('.')

        ## test OCSSW data
        if sp[0] in ['SNPP_VIIRS', 'JPSS1_VIIRS', 'JPSS2_VIIRS']:
            ## data processed from L1A from OCSSW
            ftypes = {'l1b_mod': 'L1B_MOD','geo_mod': 'GEO_MOD',
                      'l1b_img': 'L1B_IMG','geo_img': 'GEO_IMG'}
            for ft in ftypes:
                fb = [s for s in sp]
                fb[2] = ftypes[ft]

                fs = glob.glob('{}/{}'.format(dn, '.'.join(fb)))
                if len(fs) == 1:
                    if files is None: files = {}
                    files[ft] = fs[0]

        ss = sp[0][0:3]
        sl = sp[0][3:5]
        si = sp[0][5:]

        ## test LAADS DAAC data
        if ss in ['VNP', 'VJ1', 'VJ2']:
            ## data from LAADS DAAC https://search.earthdata.nasa.gov/portal/idn/search?q=
            ftypes = {'l1b_mod': '02MOD','geo_mod': '03MOD',
                      'l1b_img': '02IMG','geo_img': '03IMG'}
            for ft in ftypes:
                fs = glob.glob('{}/{}*.nc'.format(dn, '{}*{}.'.format(ss, ftypes[ft])+'.'.join(bn.split('.')[1:4])))
                if len(fs) == 1:
                    if files is None: files = {}
                    files[ft] = fs[0]

    return(files)
