## def bundle_test
## find files from AVHRR L1B/1C European Data Set bundle
## written by Quinten Vanhellemont, RBINS
## 2024-11-14

def bundle_test(bundle, extract = True):
    import os, glob, tarfile

    ## get bundle base name
    bn = os.path.basename(bundle)[0:-6]

    ## find metadata file
    meta_file = glob.glob('{}/{}.MD.XML'.format(bundle, bn))
    if len(meta_file) > 0:
        meta_file = meta_file[0]
    else:
        meta_file = None

    ## find image file
    image_file = glob.glob('{}/{}/image.nc'.format(bundle, bn))
    if len(image_file) > 0:
        image_file = image_file[0]
    else:
        image_file = None

    ## extract TAR file if not yet done
    if extract:
        tar_file = '{}/{}.TAR'.format(bundle, bn)
        if image_file is None:
            with tarfile.open(tar_file) as f:
                files = f.getnames()
                f.extractall(bundle)

            ## look for image.nc again
            image_file = glob.glob('{}/{}/image.nc'.format(bundle, '*'))
            if len(image_file) > 0:
                image_file = image_file[0]
            else:
                image_file = None

    return(image_file, meta_file)
