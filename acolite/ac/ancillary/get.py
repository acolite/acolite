## ancillary_get
## downloads and interpolates ancillary data from the ocean data server
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-10-18
## modifications: 2017-10-24 (QV) added kind keyword
##                2017-12-05 (QV) added test if data is missing
##                2018-07-18 (QV) changed acolite import name
##                2018-11-19 (QV) added verbosity option
##                2021-03-01 (QV) simplified for acg, renamed from ancillary_get

def get(date, lon, lat, local_dir=None, quiet=True, kind='linear', verbosity=0):
    import acolite as ac
    import dateutil.parser, datetime
    import os

    if type(date) == str:
        dt = dateutil.parser.parse(date)
    else:
        dt = date

    diff = (dateutil.parser.parse(datetime.datetime.now().strftime('%Y-%m-%d')) -\
            dateutil.parser.parse(dt.strftime('%Y-%m-%d'))).days
    if diff < 7:
        print('Scene too recent to get ancillary data')
        return({})

    isodate = dt.isoformat()[0:10]
    ftime = dt.hour + dt.minute/60 + dt.second/3600

    ## list and download files
    anc_local = ac.ac.ancillary.download(date = isodate, verbosity=verbosity)

    ## get ozone file
    ### use toast as fallback
    ozone_file = None
    for t in ["AURAOMI", "EPTOMS", "TOAST", 'N7TOMS']:
        for i, j in enumerate(anc_local):
            if (ozone_file is None) & (t in j):
                ozone_file = anc_local[i]

    ## set up ancillary
    anc = {'date':date, 'lon':lon, 'lat': lat, 'ftime':ftime}

    ## interpolate ozone
    if ozone_file is None:
        if verbosity > 0: print('No ozone file found for {}'.format(date))
    else:
        if os.path.exists(ozone_file):
            if verbosity > 1: print('Reading ozone from {}'.format(ozone_file))
            anc_ozone = ac.ac.ancillary.interp_ozone(ozone_file, lon, lat, kind=kind)
            for k in anc_ozone.keys(): anc[k] = anc_ozone[k]

    ## get ncep MET files
    ncep_files = [anc_local[i] for i, j in enumerate(anc_local) if "NCEPR2" in j]
    ncep_files = [f for f in ncep_files if os.path.exists(f)]

    ## interpolate MET
    if len(ncep_files) == 0:
        if verbosity > 0: print('No NCEP files found for {}'.format(date))
    else:
        if verbosity > 1: print('Reading {} ncep met files'.format(len(ncep_files)))
        anc_met = ac.ac.ancillary.interp_met(ncep_files,  lon, lat, ftime, kind=kind)
        for k in anc_met.keys(): anc[k] = anc_met[k]

    return(anc)
