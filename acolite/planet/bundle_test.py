## def bundle_test
## lists files in given directory and returns dict with band and file path
## written by Quinten Vanhellemont, RBINS
## 2018-03-12
## modifications: 2018-03-14 (QV) added option to give .tif or metadata.xml files
##                2018-03-14 (QV) improved filtering to include clipped files
##                2018-03-19 (QV) added MS files in filtering
##                2021-02-24 (QV) new version for acg
##                2022-01-11 (QV) added AnalyticMS_8b
##                2022-02-21 (QV) added Skysat, include support for unzipped API downloads

def bundle_test(bundle_in):
    import os

    if os.path.isdir(bundle_in):
        bundle = bundle_in
    else:
        bundle = os.path.dirname(bundle_in)

    ## check files in bundle
    files = []
    for f in os.listdir(bundle): files.append(os.path.join(bundle, f))
    ## include files/analytic_udm2 directory contents
    for dname in ['files', 'analytic', 'analytic_udm2', 'analytic_8b_udm2']:
        files_dir = os.path.join(bundle, dname)
        if os.path.exists(files_dir):
            for f in os.listdir(files_dir): files.append(os.path.join(files_dir, f))
    files.sort()

    datafiles = {}
    for i, file in enumerate(files):
        fname = os.path.basename(file)
        fn,ext = os.path.splitext(fname)

        if ext not in ['.json', '.tif', '.xml']: continue
        band,clp=None,''
        if 'clip' in fn:
            clp='_clip'
        if ('Analytic_metadata{}.xml'.format(clp) in fname)|\
           ('AnalyticMS_metadata{}.xml'.format(clp) in fname)|\
           ('AnalyticMS_8b_metadata{}.xml'.format(clp) in fname):
            band = 'metadata'
        if ('metadata.json' in fname):
            band = 'metadata_json'
        if ('Analytic{}.tif'.format(clp) in fname)|\
           ('AnalyticMS{}.tif'.format(clp) in fname)|\
           ('AnalyticMS_8b{}.tif'.format(clp) in fname)|\
           ('analytic{}.tif'.format(clp) in fname):
            band = 'analytic'
        if ('DN_udm{}.tif'.format(clp) in fname):
            band = 'udm'
        if ('udm2{}.tif'.format(clp) in fname):
            band = 'udm2'
        if ('analytic_dn{}.tif'.format(clp) in fname):
            band = 'analytic_dn'
        if ('panchromatic_dn{}.tif'.format(clp) in fname):
            band = 'pan_dn'
        if ('pansharpened{}.tif'.format(clp) in fname):
            band = 'pansharpened'
        if ('Analytic_SR{}.tif'.format(clp) in fname):
            band = 'sr'

        if band is None: continue
        if os.path.isfile(file):
            datafiles[band] = {"path":file, "fname":fname, "ext": ext}
    return(datafiles)
