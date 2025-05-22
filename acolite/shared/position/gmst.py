## def gmst
## computes Greenwich Mean Sidereal Time for a given time
##
## written by Quinten Vanhellemont, RBINS
## based on notes from Satellite Times by Dr. T.S. Kelso (https://celestrak.com/columns/)
## 2022-02-09
## modifications: 2025-05-22 (QV) split from eci_geometry

def gmst(time):
    import datetime, dateutil.parser
    import numpy as np
    
    ## time stuff
    if type(time) is str:
        dtime = dateutil.parser.parse(time)
    elif type(time) == datetime.datetime:
        dtime = time
    else:
        print('Provide time as string or as datetime object.')
        return

    ## find julian date
    ## du is days since JD 2451545.0 (midday UT 1st Jan 2000)
    y = dtime.year - 1
    A = int(y / 100)
    B = 2 - A + int(A / 4)
    julday_year = int(365.25 * y) + int(30.6001 * 14) + 1720994.5 + B
    doy = int(dtime.strftime('%j')) ## day of year
    du = julday_year + doy

    ## get fraction of day (out of 24 hours)
    dayfrac = dtime - datetime.datetime(dtime.year, dtime.month, dtime.day, tzinfo=dtime.tzinfo)
    deltat = dayfrac.seconds/86400

    ## get hour angle
    jdi = du + deltat
    ut = jdi + 0.5
    ut = ut - int(ut)
    jd = jdi - ut

    tu  = (jd - 2451545.0) / 36525
    gmst = 24110.54841 + tu * (8640184.812866 + tu * (0.093104 - tu * 6.2E-6))
    gmst = np.mod((gmst + 86400.0 * 1.00273790934 * ut), 86400.0)
    theta_g = np.pi * 2 * gmst / 86400.0

    return(theta_g)
