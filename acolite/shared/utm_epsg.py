## def get utm zone ant epsg for given lon lat combination
##
## written by Quinten Vanhellemont
## 2022-01-05
## modifications: 2021-02-27 (QV) tm can be datetime object


def utm_epsg(lon, lat):
    import numpy as np
    utm_zone = '{:.0f}'.format(((np.floor((float(lon) + 180) / 6 ) % 60) + 1)).zfill(2)
    if float(lat) >0:
        epsg = '326{}'.format(utm_zone)
    else:
        epsg = '327{}'.format(utm_zone)
    return(utm_zone, epsg)
