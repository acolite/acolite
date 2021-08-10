## def bundle tests
## finds image files and metadata file for DESIS bundle
## written by Quinten Vanhellemont, RBINS
## 2021-08-10

def bundle_test(bundle):
    import os, glob

    metafile = glob.glob('{}/*METADATA.xml'.format(bundle))[0]
    imagefile = glob.glob('{}/*SPECTRAL_IMAGE.tif'.format(bundle))[0]

    return(metafile, imagefile)
