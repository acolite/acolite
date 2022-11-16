## def bundle_test
## bundle test for EnMAP data
## returns dict with paths to metadata file and band files
## written by Quinten Vanhellemont, RBINS
## 2022-09-19
## modifications: 2022-09-20 (QV) changed file type handling

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

        sp = bn.split('-')
        if len(sp) == 4:
            sen, lev, scene, ftype = sp
        else:
            continue

        if ext.lower() == '.xml':
            if ftype.lower() == 'metadata':
                files_dict[ftype] = '{}'.format(fpath)
        elif ext.lower() == '.tif':
            files_dict[ftype] = '{}'.format(fpath)
    return(files_dict)
