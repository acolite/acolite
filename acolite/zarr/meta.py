## def meta
## return metadata from zarr path
## written by Quinten Vanhellemont, RBINS
## 2025-11-04
## modifications:

def meta(z):
    import os,  zarr
    md = None

    ## if str assume it is a path that can be opened
    if type(z) is str:
        f = zarr.open(z, mode='r')
        md = f.metadata.to_dict()
        f = None
    ## if zarr Group
    elif type(z) == zarr.core.group.Group:
        md = z.metadata.to_dict()
    else:
        print('Inputfile not zarr {}'.format(type(z)))

    return(md)
