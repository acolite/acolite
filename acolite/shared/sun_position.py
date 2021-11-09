## def sun position
## calculate sun position for date and time
##
## modifications: QV 2020-04-01 modified for arrays of lon and lat
##                QV 2021-11-04 return dict with more parameters

def sun_position(date, lon, lat):
    import numpy as np
    import dateutil.parser

    if type(date) is str:
        dtime = dateutil.parser.parse(date)
    else:
        dtime = date

    doy = float(dtime.strftime('%j'))

    lon = np.asarray(lon)
    if lon.shape == (): lon.shape = (1)
    lat = np.asarray(lat)
    if lat.shape == (): lat.shape = (1)

    dtor = np.pi/180
    # Get Julian date - 2400000
    ftime = dtime.hour + dtime.minute/60. + dtime.second / 3600.
    delta = dtime.year - 1949
    leap = np.floor(delta / 4.) ;# former leapyears
    jd = 32916.5 + delta * 365. + leap + doy + ftime / 24.

    # The input to the Atronomer's almanach is the difference between
    # the Julian date and JD 2451545.0 (noon, 1 January 2000)
    almanach_time = jd - 51545.

    # Ecliptic coordinates
    # Mean longitude
    mnlong = 280.460 + .9856474 * almanach_time
    mnlong = mnlong % 360.
    if mnlong<0:mnlong+=360.

    # Mean anomaly
    mnanom = 357.528 + .9856003 * almanach_time
    mnanom = mnanom % 360.
    if mnanom<0:mnanom+=360.
    mnanom = mnanom * dtor

    distance = 1.00014 - 0.01671 * np.cos(mnanom) - 0.00014 * np.cos(2.*mnanom)

    # Ecliptic longitude and obliquity of ecliptic
    eclong = mnlong + 1.915 * np.sin(mnanom) + 0.020 * np.sin(2 * mnanom)
    eclong = eclong % 360.
    if eclong<0:eclong+=360.
    oblqec = 23.429 - 0.0000004 * almanach_time
    oblqec = 23.439 - 0.0000004 * almanach_time
    eclong = eclong * dtor
    oblqec = oblqec * dtor

    # Celestial coordinates
    # Right ascension and declination
    num = np.cos(oblqec) * np.sin(eclong)
    den = np.cos(eclong)
    ra = np.arctan(num / den)
    if den < 0: ra+=np.pi
    if (den > 0) & (num > 0): ra+= np.pi * 2.
    dec = np.arcsin(np.sin(oblqec) * np.sin(eclong))

    # Local coordinates
    # Greenwich mean sidereal time
    gmst = 6.697375 + .0657098242 * almanach_time + ftime
    gmst = gmst % 24.
    if gmst < 0: gmst+=24.

    # Local mean sidereal time
    lmst = (lon / 15.) + np.asarray(gmst)
    lmst = lmst % 24.
    lmst[lmst<0]+=24.
    lmst *= 15. * dtor

    # Hour angle
    ha = lmst - ra
    ha[ha < -1*np.pi] += np.pi * 2.
    ha[ha > np.pi] -= np.pi * 2.

    # Azimuth and elevation
    el = np.arcsin(np.sin(dec) * np.sin(lat*dtor) + np.cos(dec) * np.cos(lat*dtor) * np.cos(ha))
    az = np.arcsin(-1*np.cos(dec) * np.sin(ha) / np.cos(el))

    cosAzPos = 0 < (np.sin(dec) - np.sin(el) * np.sin(lat*dtor))
    sinAzNeg = np.sin(az) < 0
    az[cosAzPos & sinAzNeg] += np.pi * 2
    az[cosAzPos == False] = np.pi - az[cosAzPos == False]
    az[np.isinf(az)] = np.pi/2

    az/=dtor
    zen = 90.-el/dtor

    d = {}
    d['julian_date'] = jd
    d['almanach_time'] = almanach_time
    d['mean_longitude'] = mnlong
    d['mean_anomaly'] = mnanom
    d['elevation'] = el/dtor
    d['zenith'] = zen
    d['azimuth'] = az
    d['distance'] = distance
    d['right_ascension'] = ra
    d['declination'] = dec
    d['greenwich_mean_sidereal_time'] = gmst
    d['local_mean_sidereal_time'] = lmst
    d['hour_angle'] = ha
    return(d)
