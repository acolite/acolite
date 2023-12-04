## def bundle_test
## bundle test for Haiyang data
## returns dict with paths to metadata file and band files
## written by Quinten Vanhellemont, RBINS
## 2023-02-19
## modifications: 2023-12-04 (QV) added HDF support

def bundle_test(bundle):
    import os, glob
    if os.path.isdir(bundle):
        files = glob.glob('{}/*'.format(bundle))
    else:
        files = glob.glob('{}/*'.format(os.path.dirname(bundle)))

    files.sort()
    files_dict = {}
    for file in files:
        fpath = os.path.abspath(file)
        bn = os.path.basename(file)
        bn, ext = os.path.splitext(bn)
        if ext.lower() == '.tiff':
            ftype = 'image'
        elif ext.lower() == '.xml':
            ftype = 'metadata'
        elif ext.lower() == '.h5':
            ftype = 'image'
        else: continue
        files_dict[ftype] = '{}'.format(fpath)
    return(files_dict)
