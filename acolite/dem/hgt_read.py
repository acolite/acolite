## def hgt_read
## reads DEM HGT SRTM files
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-07-17
##                2019-04-24 (QV) added support for zip files
##                2021-04-07 (QV) changed numpy import

def hgt_read(file):
    import struct
    import numpy as np

    if '.gz' in file:
        import gzip
        with gzip.open(file,'rb') as f:
            data_read = f.read()
    elif '.zip' in file:
        import zipfile, os
        zfile = '{}.{}'.format(os.path.basename(file).split('.')[0], 'hgt')
        with zipfile.ZipFile(file, mode='r') as f:
            data_read = f.read(zfile)
    else:
        with open(file,'rb') as f:
            data_read = f.read()

    dim = (1201,1201)

    ## big endian, unsigned shorts
    data = struct.unpack('>{}'.format('H'*dim[0]*dim[1]),data_read)
    data = np.asarray(data).reshape(dim)

    data[data > 32768] -= 65535
    return(data)
