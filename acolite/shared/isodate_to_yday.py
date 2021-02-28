## def isodate_to_yday
## converts isodate (YYYY-MM-DDTHH:mm:ss) to day of year (fractional)
## written by Quinten Vanhellemont
## 2018-05-15
## modifications: 2021-02-27 (QV) tm can be datetime object

def isodate_to_yday(tm, return_yf=False):
    import dateutil.parser
    from datetime import datetime

    if type(tm) is str:
        tm = dateutil.parser.parse(tm)

    y = tm.year

    ## get day fraction
    hour = tm.hour + tm.minute/60 + tm.second/3600
    df = hour/24.

    ## get year fraction
    ylen = ((datetime(y+1, 1, 1) - datetime(y, 1, 1)).days)
    doy = float(tm.strftime('%j'))
    doy+=df

    yf = doy/ylen

    if return_yf:
        return(y+yf)
    else:
        return(tm, y, yf)
