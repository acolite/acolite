## def bundle_test
## lists files in given directory and returns dict with band and file path
## written by Quinten Vanhellemont, RBINS
## 2018-03-12
## modifications: 2018-03-14 (QV) added option to give .tif or metadata.xml files
##                2018-03-14 (QV) improved filtering to include clipped files
##                2018-03-19 (QV) added MS files in filtering
##                2021-02-24 (QV) new version for acg
##                2022-01-11 (QV) added AnalyticMS_8b

def bundle_test(bundle_in):
    import os

    if os.path.isdir(bundle_in):
        bundle = bundle_in
    else:
        bundle = os.path.dirname(bundle_in)

    files = os.listdir(bundle)
    datafiles = {}
    for i, fname in enumerate(files):
        fn,ext = os.path.splitext(fname)
        if ext == '.json': continue
        if ext not in ['.tif', '.xml']: continue
        band,clp=None,''
        if 'clip' in fn:
            clp='_clip'
        if ('Analytic_metadata{}.xml'.format(clp) in fname)|\
           ('AnalyticMS_metadata{}.xml'.format(clp) in fname)|\
           ('AnalyticMS_8b_metadata{}.xml'.format(clp) in fname):
            band = 'metadata'
        if ('Analytic{}.tif'.format(clp) in fname)|\
           ('AnalyticMS{}.tif'.format(clp) in fname)|\
           ('AnalyticMS_8b{}.tif'.format(clp) in fname):
            band = 'analytic'
        if ('DN_udm{}.tif'.format(clp) in fname):
            band = 'udm'
        if ('Analytic_SR{}.tif'.format(clp) in fname):
            band = 'sr'

        if band is None: continue
        file = '{}/{}'.format(bundle,fname)
        if os.path.isfile(file):
            datafiles[band] = {"path":file, "fname":fname}
    return(datafiles)
