## def multi_tile_extent
## gets dct for multiple tiles
##
## written by Quinten Vanhellemont, RBINS
## 2024-03-28
## modifications:

def multi_tile_extent(inputfile, dct = None):
    import acolite as ac
    import numpy as np
    import glob

    if type(inputfile) is not list:
        inputfile_ = [inputfile]
    else:
        inputfile_ = [f for f in inputfile]

    ## get the projection dcts for all scenes
    dct_list = []
    for bundle in inputfile_:
        ## get scene projection and extent
        mtl = glob.glob('{}/{}'.format(bundle, '*MTL.txt'))
        ## add ALI MTL files
        mtl += glob.glob('{}/{}'.format(bundle, '*MTL_L1T.TXT'))
        mtl += glob.glob('{}/{}'.format(bundle, '*MTL_L1GST.TXT'))

        if len(mtl) == 0: continue
        meta = ac.landsat.metadata_read(mtl[0])
        dct_list.append(ac.landsat.projection(meta))
    ## merge the dcts
    dct = ac.shared.projection_merge(dct_list, dct = dct)

    ## add one to fix Landsat dimensions
    dct['xdim']+=1
    dct['ydim']+=1
    dct['dimensions'] = dct['xdim'], dct['ydim']

    return(dct)
