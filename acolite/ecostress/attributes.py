## def attributes
## reads attributes from ECOSTRESS HDF file
## written by Quinten Vanhellemont, RBINS
## 2022-08-11
## modifications: 2022-08-16 (QV) convert to str

def attributes(file):
    import h5py
    h5_gatts = {}
    with h5py.File(file, mode='r') as f:
        for a in f['StandardMetadata'].keys():
            h5_gatts[a] = f['StandardMetadata'][a][()]
        for a in f['L1GEOMetadata'].keys():
            h5_gatts[a] = f['L1GEOMetadata'][a][()]

    for a in h5_gatts:
        if type(h5_gatts[a]) == bytes:
            h5_gatts[a] = h5_gatts[a].decode('utf-8')

    return(h5_gatts)
