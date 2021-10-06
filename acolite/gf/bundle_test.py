## def bundle tests
## finds image files and metadata file for GF bundle
## written by Quinten Vanhellemont, RBINS
## 2021-08-09

def bundle_test(bundle):
    import os, glob
    files = glob.glob('{}/*'.format(bundle))
    files.sort()

    info = {}
    for f in files:
        bn,  ex = os.path.splitext((os.path.basename(f)))
        if ('_thumb' in bn) or ('.rpb.aux' in bn): continue
        if ('order' == bn) & (ex == '.xml'): continue
        if bn not in info: info[bn] = {}
        info[bn][ex] = {'path': f}

    tiles = []
    metafile = None
    for k in info:
        #print(k)
        keys = list(info[k].keys())
        if ('.xml' in keys) & ('.tiff' not in keys):
            metafile = info[k]['.xml']['path']
        if ('.tiff' in keys):
            tiles.append(info[k]['.tiff']['path'])

    return(tiles, metafile)
