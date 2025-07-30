## def bundle_test
## finds Huanjing files
##
## written by Quinten Vanhellemont, RBINS
## 2025-07-29
## modifications:

def bundle_test(bundle):
    import glob, os

    if os.path.isdir(bundle):
        bundle_files = glob.glob('{}/HJ2*'.format(bundle))
    else:
        bundle_files = glob.glob('{}/HJ2*'.format(os.path.dirname(bundle)))
    bundle_files.sort()

    ## find tiff and xml
    imagefile, metafile = None, None
    for file in bundle_files:
        bn = os.path.basename(file)
        dn = os.path.dirname(file)

        bn, ext = os.path.splitext(bn)

        if ext == '.xml':
            metafile = '{}'.format(file)
        if ext == '.tiff':
            imagefile = '{}'.format(file)
    return(imagefile, metafile)
