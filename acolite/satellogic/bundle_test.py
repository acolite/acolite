## def bundle_test
## identifies files for Satellogic NewSat data
## note: works with L1C sample data from s3://satellogic-sample-data
## written by Quinten Vanhellemont, RBINS
## 2026-11-17
## modifications: 2026-11-18 (QV) added support for passing a scene file directly

def bundle_test(bundle):
    import os, glob

    ## get image and metadata file
    if os.path.isdir(bundle):
        files = glob.glob('{}/*'.format(bundle))
    else:
        dn = os.path.dirname(bundle)
        bn = os.path.basename(bundle)
        files = glob.glob('{}/{}*'.format(dn, bn[0:bn.find('_SN')]))
    files.sort()
    basenames = [os.path.basename(f) for f in files]

    items = {}
    for f in files:
        bn = os.path.basename(f)
        bn = bn[0:bn.find('.')]
        it = bn[0:32]
        dt = bn[32:]
        if it not in items: items[it] = {}
        items[it][dt] = f

    ## check item presence
    delete = []
    for it in items:
        if ('analytic' not in items[it]) | ('2_2_0' not in items[it]):
            delete.append(it)
    for it in delete: del items[it]

    if len(items) == 0: return
    return(items)
