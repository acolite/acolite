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
##                2022-11-17 (QV) delete series by default
##                2023-10-16 (QV) added GMAO data, rescale and rename parameters in this function
##                2023-12-28 (QV) added GMAO_IT

def get(date, lon, lat, local_dir=None, quiet=True, kind='linear', verbosity=0, keep_series=False):
    import acolite as ac
    import dateutil.parser, datetime
    import os

    if type(date) == str:
        dt = dateutil.parser.parse(date)
    else:
        dt = date

    isodate = dt.isoformat()#[0:10]
    ftime = dt.hour + dt.minute/60 + dt.second/3600

    diff = (dateutil.parser.parse(datetime.datetime.now().strftime('%Y-%m-%d')) -\
            dateutil.parser.parse(dt.strftime('%Y-%m-%d'))).days
    # if diff < 7:
    #     print('Scene too recent to get ancillary data: {}'.format(isodate))
    #     return({})

    if isodate < '1978-10-27':
        print('Scene too old to get ancillary data: {}'.format(isodate))
        return({})

    ## list and download files
    anc_local = ac.ac.ancillary.download(date = isodate, verbosity=verbosity)

    ## find if we have merra2 files
    gmao_files = [file for file in anc_local if ('GMAO_MERRA2' in os.path.basename(file)) & (os.path.exists(file))]
    if len(gmao_files) == 0: gmao_files = [file for file in anc_local if ('GMAO_FP' in os.path.basename(file)) & (os.path.exists(file))]
    if len(gmao_files) == 0: gmao_files = [file for file in anc_local if ('GMAO_IT' in os.path.basename(file)) & (os.path.exists(file))]

    if len(gmao_files) == 2:
        if verbosity > 1:
            print('Using GMAO GEOS ancillary data:')
            for file in gmao_files: print(file)

        ## set up ancillary
        anc = {'date':date, 'lon':lon, 'lat': lat, 'ftime':ftime, 'type': 'merra2', 'data': {}}
        anc_gmao = ac.ac.ancillary.interp_gmao(gmao_files,  lon, lat, isodate, method=kind)
        for k in anc_gmao.keys():
            if (not keep_series) & ('series' in anc_gmao[k]): del anc_gmao[k]['series']
            anc['data'][k] = anc_gmao[k]

    else:
        ## get ozone file
        ### use toast as fallback
        ozone_file = None
        for t in ["AURAOMI", "EPTOMS", "TOAST", 'N7TOMS']:
            for i, j in enumerate(anc_local):
                if (ozone_file is None) & (t in j) & (os.path.exists(anc_local[i])):
                    ozone_file = anc_local[i]

        ## set up ancillary
        anc = {'date':date, 'lon':lon, 'lat': lat, 'ftime':ftime, 'data': {}}

        ## interpolate ozone
        if ozone_file is None:
            if verbosity > 0: print('No ozone file found for {}'.format(date))
        else:
            if os.path.exists(ozone_file):
                if verbosity > 1: print('Reading ozone from {}'.format(ozone_file))
                anc_ozone = ac.ac.ancillary.interp_ozone(ozone_file, lon, lat, kind=kind)
                for k in anc_ozone.keys():
                    if (not keep_series) & ('series' in anc_ozone[k]): del anc_ozone[k]['series']
                    anc['data'][k] = anc_ozone[k]

        ## get ncep MET files
        ncep_files = [anc_local[i] for i, j in enumerate(anc_local) if "NCEPR2" in j]
        if len(ncep_files) == 0:
            ncep_files = [anc_local[i] for i, j in enumerate(anc_local) if "NCEP" in j]
        ncep_files = [f for f in ncep_files if os.path.exists(f)]

        ## interpolate MET
        if len(ncep_files) == 0:
            if verbosity > 0: print('No NCEP files found for {}'.format(date))
        else:
            if verbosity > 1: print('Reading {} ncep met files'.format(len(ncep_files)))
            anc_met = ac.ac.ancillary.interp_met(ncep_files,  lon, lat, ftime, kind=kind)
            for k in anc_met.keys():
                if (not keep_series) & ('series' in anc_met[k]): del anc_met[k]['series']
                anc['data'][k] = anc_met[k]
        ## add type
        anc['type'] = 'ncep'

    ## rescale ancillary data
    anc_name = ['uoz', 'uwv', 'z_wind', 'm_wind', 'pressure']
    if anc['type'] == 'ncep':
        anc_keys = ['ozone', 'p_water', 'z_wind', 'm_wind', 'press']
        anc_fact = [1./1000., 1./10., 1., 1., 1.]
    elif anc['type'] == 'merra2':
        anc_keys = ['TO3', 'TQV', 'U10M', 'V10M', 'PS']
        anc_fact = [1./1000., 1./10., 1., 1., 1./100.]
    for i, k in enumerate(anc_keys):
        if k not in anc['data']: continue
        anc[anc_name[i]] = anc['data'][k]['interp'] * anc_fact[i]
    ## compute wind speed
    if ('z_wind' in anc) & ('m_wind' in anc):
        anc['wind'] = (((anc['z_wind'])**2) + ((anc['m_wind'])**2))**0.5

    return(anc)
