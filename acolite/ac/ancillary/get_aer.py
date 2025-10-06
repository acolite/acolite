## get_aer
## downloads and interpolates ancillary aerosol data from the ocean data server
##
## written by Quinten Vanhellemont, RBINS
## 2025-10-06
## modifications: 2025-10-06 (QV) determine file_types based on date

def get_aer(date, lon, lat, local_dir = None,
            file_types = None, datasets = ['TOTEXTTAU', 'TOTSCATAU', 'TOTANGSTR'],
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
        diff = (dateutil.parser.parse(datetime.datetime.now().strftime('%Y-%m-%d')) -\
                dateutil.parser.parse(dt.strftime('%Y-%m-%d'))).days
        file_types = ['GMAO_MERRA2_AER']
        if diff < 40: file_types = ['GMAO_IT_AER']
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
