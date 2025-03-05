## def bundle_test
## finds xml file in WV bundle
## written by Quinten Vanhellemont, RBINS
## 2025-03-01
## modifications: 2025-03-05 (QV) commented printout

def bundle_test(bundle):
    import glob, os

    if os.path.isdir(bundle):
        search_path = '{}/'.format(bundle)
    else:
        bn, ext = os.path.splitext(bundle)
        search_path = '{}'.format(bn)

    ## search for xml files
    metafiles = glob.glob('{}{}'.format(search_path,'*.XML'))
    metafiles += glob.glob('{}{}'.format(search_path,'*.xml'))
    metafiles.sort()

    ## find non aux or readme files
    ## sorting should get MS file over PAN file
    if len(metafiles)>0:
        idx = 0
        if len(metafiles) >= 1:
            for idx, mf in enumerate(metafiles):
                if ('.aux.' not in mf) & ('README' not in mf) & ('(1)' not in mf):
                    break
        metafile = metafiles[idx]
        return(metafile)
    else:
        #print('No metadata found for {}'.format(bundle))
        return
