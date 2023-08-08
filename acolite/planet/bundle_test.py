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
##                2022-10-26 (QV) added scene_id to datafiles
##                2023-04-17 (QV) fix for Skysat scene_ids and selecting one from multiple files
##                2023-04-18 (QV) added PSScene and files/PSScene dname options
##                                added support for NTF files
##                2023-05-08 (QV) added support for composite files
##                2023-05-25 (QV) set sid to None if manifest is given

def bundle_test(bundle_in):
    import os, glob

    if os.path.isdir(bundle_in):
        bundle = bundle_in
        sid = None
    else:
        bundle = os.path.dirname(bundle_in) #+ os.path.sep
        fn = os.path.basename(bundle_in)
        if fn == 'manifest.json':
            sid = None
        else:
            sid = fn[0:23]
            if 'ssc' in fn: sid = fn[0:27]

    ## check files in bundle
    files = []
    for f in os.listdir(bundle): files.append(os.path.join(bundle, f))
    ## include files/analytic_udm2 directory contents
    for dname in ['files', 'analytic', 'analytic_udm2', 'analytic_8b_udm2', \
                  'PSScene', 'files/PSScene']:
        files_dir = os.path.join(bundle, dname)
        if os.path.exists(files_dir):
            if 'PSScene' in dname:
                for f in os.listdir(files_dir): files.append(os.path.join(files_dir, f)) ## newer zip bundles
                for f in glob.glob('{}/*/analytic*/*'.format(files_dir)): files.append(f) ## older zip bundles
            else:
                for f in os.listdir(files_dir): files.append(os.path.join(files_dir, f))
    files.sort()

    datafiles = {}
    for i, file in enumerate(files):
        fname = os.path.basename(file)
        fn,ext = os.path.splitext(fname)
        if len(fn) < 23: continue
        if sid is None:
            scene_id = fn[0:23]
            if 'ssc' in fn: scene_id = fn[0:27]
        else:
            scene_id = sid
        if scene_id not in fn: continue
        if ext not in ['.json', '.tif', '.xml', '.ntf']: continue
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
        if ('Analytic{}_file_format.ntf'.format(clp) in fname)|\
           ('AnalyticMS{}_file_format.ntf'.format(clp) in fname)|\
           ('AnalyticMS_8b{}_file_format.ntf'.format(clp) in fname)|\
           ('analytic{}_file_format.ntf'.format(clp) in fname):
            band = 'analytic_ntf'
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
        if ('Analytic_SR{}.tif'.format(clp) in fname)|\
           ('AnalyticMS_SR_8b{}.tif'.format(clp) in fname):
            band = 'sr'

        if ('composite.tif' in fname): band = 'composite'
        if ('composite_udm2.tif' in fname): band = 'composite_udm2'

        if band is None: continue
        if os.path.isfile(file):
            if scene_id not in datafiles: datafiles[scene_id] = {}
            datafiles[scene_id][band] = {"path":file, "fname":fname, "ext": ext}
    return(datafiles)
