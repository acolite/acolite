## def bundle_test
## identifies files for Satellogic NewSat data
## note: works with L1C sample data from s3://satellogic-sample-data
## written by Quinten Vanhellemont, RBINS
## 2026-11-17
## modifications: 2026-11-18 (QV) added support for passing a scene file directly
##                2026-11-19 (QV) added support for L1D files

def bundle_test(bundle):
    import os, glob

    ## get image and metadata file
    if os.path.isdir(bundle):
        search_dir = '{}'.format(bundle)
    else:
        dn = os.path.dirname(bundle)
        bn = os.path.basename(bundle)
        search_dir = '{}/{}'.format(dn, bn[0:bn.find('_SN')])

    print(search_dir)
    ## find files and rasters
    files = glob.glob('{}/*'.format(search_dir))
    files += glob.glob('{}/rasters/*'.format(search_dir))
    files.sort()

    ## identify files
    items = {}
    for f in files:
        bn = os.path.basename(f)
        if bn.endswith('.ovr'): continue
        if bn == 'rasters': continue

        bn = bn[0:bn.find('.')]
        it = bn[0:32]
        dt = bn[32:]
        if dt in ['2_2_0', 'metadata_stac']: dt = 'meta'
        if it not in items: items[it] = {}
        items[it][dt] = f

    ## check item presence
    delete = []
    for it in items:
        if ('analytic' in items[it]) & ('meta' in items[it]): continue # L1C
        if ('TOA_0' in items[it]) & ('meta' in items[it]): continue # L1D
        delete.append(it)
    for it in delete: del items[it]

    if len(items) == 0: return
    return(items)
