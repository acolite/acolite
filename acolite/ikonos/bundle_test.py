## def bundle_test
## bundle test for (old?) IKONOS data
## returns dict with paths to metadata file and band files
## written by Quinten Vanhellemont, RBINS
## 2022-09-14
## modifications: 2022-09-15 (QV) made band files a list, in case there's more components
##                2022-11-21 (QV) skip EUSI_logo

def bundle_test(bundle):
    import os, glob
    files = glob.glob('{}/*'.format(bundle))
    files.sort()
    files_dict = {}
    for file in files:
        fpath = os.path.abspath(file)
        bn = os.path.basename(file)
        bn, ext = os.path.splitext(bn)
        if bn in ['EUSI_logo']: continue

        if ext == '.txt':
            if '_metadata' in bn:
                files_dict['metadata'] = '{}'.format(fpath)
        elif ext == '.jpg':
            band = bn.split('_')[-3]
            files_dict[band] = '{}'.format(fpath)
        elif ext == '.tif':
            band = bn.split('_')[-2]
            if band not in files_dict: files_dict[band] = []
            files_dict[band].append('{}'.format(fpath))
    return(files_dict)
