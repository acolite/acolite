## def datascl
## linearly rescales data from the range [dmin,dmax] to [tmin,tmax], defaults to scaling to byte range and type
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-06-29
## modifications: 2021-02-24 (QV) small updates for acg
##

def datascl(data, dmin=None, dmax=None, tmin=0, tmax=255, dtype='uint8', percentiles=(1,99)):
    import numpy as np
    if dmin == None:
        if percentiles is not None:
            dmin = np.nanpercentile(data, percentiles[0])
        else:
            dmin = np.nanmin(data)
    if dmax == None:
        if percentiles is not None:
            dmax = np.nanpercentile(data, percentiles[1])
        else:
            dmax = np.nanmax(data)
    data=np.interp(data, [dmin,dmax],[tmin,tmax])
    if dtype != None: data=data.astype(dtype)
    return(data)
