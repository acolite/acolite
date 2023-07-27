## def bundle test
## finds image files and metadata files for SDGSAT1 KX10 bundle
## written by Quinten Vanhellemont, RBINS
## 2023-01-03
## modifications: 2023-02-18 (QV) track all tiffs

def bundle_test(bundle, sen='KX10_MII'):
    import glob, os

    ##
    if os.path.isfile(bundle) & ('.meta.xml' in bundle):
        metafiles = [bundle]
    else:
        metafiles = glob.glob('{}/{}_*.meta.xml'.format(bundle, sen))
        metafiles.sort()

    ## find associated cal and image files
    calfiles = []
    imgfiles = []
    for mf in metafiles:
        dn = os.path.dirname(mf)
        bn = os.path.basename(mf)
        sn = bn.replace('.meta.xml', '')
        im = glob.glob('{}/{}_*.tif'.format(dn, sn))#[0]
        im.sort()
        imgfiles.append(im)

        cf = '{}/{}'.format(dn, bn.replace('.meta.xml', '.calib.xml'))
        calfiles.append(cf)
    return(metafiles, calfiles, imgfiles)
