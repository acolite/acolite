## def attributes
## reads attributes from ECOSTRESS HDF file
## written by Quinten Vanhellemont, RBINS
## 2022-08-11

def attributes(file):
    import h5py
    h5_gatts = {}
    with h5py.File(file, mode='r') as f:
        for a in f['StandardMetadata'].keys():
            h5_gatts[a] = f['StandardMetadata'][a][()]
        for a in f['L1GEOMetadata'].keys():
            h5_gatts[a] = f['L1GEOMetadata'][a][()]
    return(h5_gatts)
