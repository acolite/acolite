## def zarr.multi_tile_extent
## gets dct for multiple tiles from zarr paths
##
## written by Quinten Vanhellemont, RBINS
## 2025-11-12
## modifications:

def multi_tile_extent(inputfile, dct = None, s2_target_res = 10):
    import acolite as ac

    if type(inputfile) is not list:
        inputfile_ = [inputfile]
    else:
        inputfile_ = [f for f in inputfile]

    ## get the projection dcts for all scenes
    dct_list = []
    for bundle in inputfile_:
        dct_list.append(ac.sentinel2.zarr.projection(bundle, s2_target_res = s2_target_res))

    ## merge the dcts
    dct = ac.shared.projection_merge(dct_list, dct = dct)

    return(dct)
