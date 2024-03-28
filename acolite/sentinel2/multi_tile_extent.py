## def multi_tile_extent
## gets dct for multiple tiles
##
## written by Quinten Vanhellemont, RBINS
## 2024-03-27
## modifications: 2024-03-28 (QV) moved merging to separate function

def multi_tile_extent(inputfile, dct = None):
    import acolite as ac
    import numpy as np

    if type(inputfile) is not list:
        inputfile_ = [inputfile]
    else:
        inputfile_ = [f for f in inputfile]

    ## get the projection dcts for all scenes
    dct_list = []
    for bundle in inputfile_:
        safe_files = ac.sentinel2.safe_test(bundle)
        granule = safe_files['granules'][0]
        grmeta = ac.sentinel2.metadata_granule(safe_files[granule]['metadata']['path'])
        dct_list.append(ac.sentinel2.projection(grmeta, s2_target_res=int(ac.settings['run']['s2_target_res'])))

    ## merge the dcts
    dct = ac.shared.projection_merge(dct_list, dct = dct)

    return(dct)
