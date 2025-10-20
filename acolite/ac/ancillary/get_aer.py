## get_aer
## downloads and interpolates ancillary aerosol data from the ocean data server
##
## written by Quinten Vanhellemont, RBINS
## 2025-10-06
## modifications: 2025-10-06 (QV) determine file_types based on date
##                2025-10-15 (QV) added nrt_days as keyword
##                2025-10-20 (QV) added ancillary_aerosol, ancillary_aerosol_nrt and ancillary_aerosol_nrt_days as settings


def get_aer(date, lon, lat, local_dir = None,
            nrt_days = None, file_types = None, datasets = ['TOTEXTTAU', 'TOTSCATAU', 'TOTANGSTR'],
            quiet = True, kind = 'linear', verbosity = 0, keep_series = False):

    import acolite as ac
    import dateutil.parser, datetime, os, netCDF4
    import scipy.interpolate
    import numpy as np

    ## available datasets
    if datasets is None:
        datasets = ['BCEXTTAU', 'BCSCATAU', 'DUEXTTAU', 'DUSCATAU', 'SSEXTTAU', 'SSSCATAU',
                    'SUEXTTAU', 'SUSCATAU', 'OCEXTTAU', 'OCSCATAU', 'TOTEXTTAU', 'TOTSCATAU', 'TOTANGSTR']

    if type(date) == str:
        dt = dateutil.parser.parse(date)
    elif type(date) == datetime.datetime:
        dt = date
    else:
        print('Provide date as string or as datetime object, currently: {}'.format(type(date)))
        print('date = {}'.format(date))
        return

    ## get date
    isodate = dt.isoformat()
    ftime = dt.hour + dt.minute/60 + dt.second/3600
    if isodate < '1978-10-27':
        print('Scene too old to get ancillary data: {}'.format(isodate))
        return()

    if file_types is None:
        file_types = ac.settings['run']['ancillary_aerosol']
        if type(file_types) is not list: file_types = [file_types]
        file_types_nrt = ac.settings['run']['ancillary_aerosol_nrt']
        if type(file_types_nrt) is not list: file_types_nrt = [file_types_nrt]

        ## find out if NRT data are to be used
        nrt_days = ac.settings['run']['ancillary_aerosol_nrt_days']
        diff = (dateutil.parser.parse(datetime.datetime.now().strftime('%Y-%m-%d')) -\
                dateutil.parser.parse(dt.strftime('%Y-%m-%d'))).days
        if diff < nrt_days: file_types = file_types_nrt
        print('Using file types {} for date {} days before today.'.format(file_types, diff))

    ## list and download files
    anc_local = ac.ac.ancillary.download(date = isodate, verbosity = verbosity, file_types = file_types)

    ## set up dict
    anc = {'date':date, 'lon':lon, 'lat': lat, 'ftime':ftime, 'type': file_types, 'data': {}}

    ## interpolate to lon/lat and time
    if len(anc_local) > 0:
        anc_gmao = ac.ac.ancillary.interp_gmao(anc_local,  lon, lat, isodate, method = kind, datasets = datasets)
        for k in anc_gmao.keys():
            if (not keep_series) & ('series' in anc_gmao[k]): del anc_gmao[k]['series']
            anc['data'][k] = anc_gmao[k]

    return(anc)
