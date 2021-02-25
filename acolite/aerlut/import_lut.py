## def import_lut
## imports LUT made with 6SV and converts to NetCDF
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07-05
## modifications:   2020-07-14 (QV)
##                  2021-01-16 (QV) added support for bz2 compressed luts
##                  2021-02-24 (QV) removed obsolete code
##                  2021-02-25 (QV) changed position of lut files (removed lutid directory), added removal of unzipped file

def import_lut(lutid,lutdir,override=0):
    import os, sys
    import numpy as np

    lutnc=lutdir+'/'+lutid+'.nc'
    lut = None

    ## extract bz2 files
    unzipped = False
    lutncbz2 = '{}.bz2'.format(lutnc)
    if (not os.path.isfile(lutnc)) & (os.path.isfile(lutncbz2)):
        import bz2, shutil
        unzipped = True
        with bz2.BZ2File(lutncbz2) as fi, open(lutnc,"wb") as fo:
            shutil.copyfileobj(fi,fo)
    ## end extract bz2 files

    ## read dataset from NetCDF
    try:
        if os.path.isfile(lutnc):
            from netCDF4 import Dataset
            nc = Dataset(lutnc)
            meta=dict()
            for attr in nc.ncattrs():
                attdata = getattr(nc,attr)
                if isinstance(attdata,str): attdata = attdata.split(',')
                meta[attr]=attdata
            lut = nc.variables['lut'][:]
            nc.close()
    except:
        print(sys.exc_info()[0])
        print('Failed to open LUT data from NetCDF (id='+lutid+')')

    if unzipped: os.remove(lutnc) ## clear unzipped LUT

    if lut is None:
        print('Could not import LUT {} from {}'.format(lutid, lutdir))
        return()

    ## for the  and Continental and Urban models (1,3)
    ## romix nans were retrieved for wavelengths > 2 micron and aot == 0.001
    ## 500mb for C+U and 1013/1100 for U
    ## if any nans set then to 0
    sub = np.where(np.isnan(lut))
    lut[sub] = 0

    return(lut, meta)
