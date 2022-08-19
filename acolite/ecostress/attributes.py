## def attributes
## reads attributes from ECOSTRESS HDF file
## written by Quinten Vanhellemont, RBINS
## 2022-08-11
## modifications: 2022-08-16 (QV) convert to str
##                2022-08-19 (QV) test if att exist

def attributes(file):
    import h5py
    h5_gatts = {}
    with h5py.File(file, mode='r') as f:
        for k in ['StandardMetadata', 'L1B_RADMetadata', 'L1GEOMetadata']:
            if k in f:
                for a in f[k].keys():
                    h5_gatts[a] = f[k][a][()]
    for a in h5_gatts:
        if type(h5_gatts[a]) == bytes:
            h5_gatts[a] = h5_gatts[a].decode('utf-8')

    return(h5_gatts)
