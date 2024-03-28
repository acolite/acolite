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
    for ifile in inputfile_:
        bundle, files, fk = ifile

        image_file = None
        if 'analytic' in files:
            image_file = files['analytic']['path']
        elif 'analytic_ntf' in files:
            image_file = files['analytic_ntf']['path']
        elif 'pansharpened' in files:
            image_file = files['pansharpened']['path']
        elif 'composite' in files:
            image_file = files['composite']['path']

        sr_image_file = None
        if 'sr' in files: sr_image_file = files['sr']['path']

        if image_file is not None:
            dct_ = ac.shared.projection_read(image_file)
        elif sr_image_file is not None:
            dct_ = ac.shared.projection_read(sr_image_file)
        else:
            continue

        dct_list.append(dct_)

    ## merge the dcts
    dct = ac.shared.projection_merge(dct_list, dct = dct)

    ## add one to fix Landsat dimensions
    #dct['xdim']+=1
    #dct['ydim']+=1
    #dct['dimensions'] = dct['xdim'], dct['ydim']

    return(dct)
