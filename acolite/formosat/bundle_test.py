## def bundle test
## finds image files and metadata file for FORMOSAT bundle
## written by Quinten Vanhellemont, RBINS
## 2022-04-12

def bundle_test(bundle):
    import os, glob
    files = glob.glob('{}/*'.format(bundle))
    files.sort()

    info = {}
    for f in files:
        bn,  ex = os.path.splitext((os.path.basename(f)))
        if bn not in info: info[bn] = {}
        info[bn][ex] = {'path': f}

    tiles = []
    metafile = None
    for k in info:
        #print(k)
        keys = list(info[k].keys())
        if ('.dim' in keys):
            metafile = info[k]['.dim']['path']
        if ('.tif' in keys):
            tiles.append(info[k]['.tif']['path'])

    return(tiles, metafile)
