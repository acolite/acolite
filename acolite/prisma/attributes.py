## Read PRISMA HDF attributes
## QV 2021-07-14

def attributes(file):
    import h5py
    h5_gatts = {}
    with h5py.File(file, mode='r') as f:
        for a in f.attrs.keys(): h5_gatts[a] = f.attrs[a]
    return(h5_gatts)
