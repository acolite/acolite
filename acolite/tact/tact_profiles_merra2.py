# def tact_profiles_merra2
## gets and reformats merra2 profiles for a given limit and dataset
## currently uses ASCII interface which is not ideal
##
## more info on this dataset at https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/
## MERRA-2 inst3_3d_asm_Np: 3d,3-Hourly,Instantaneous,Pressure-Level,Assimilation,Assimilated Meteorological Fields V5.12.4
##
## written by Quinten Vanhellemont, RBINS
## 2023-02-02
## modifications: 2023-08-07 (QV) moved url_base to ACOLITE config file

def tact_profiles_merra2(isotime, limit, obase = None, override = False, verbosity = 5, url_base = None):
    print(isotime, limit)

    import os, json
    import time, dateutil.parser, datetime
    import requests, netrc
    import numpy as np

    import acolite as ac

    ## load auth
    nr = netrc.netrc()
    ret = nr.authenticators('earthdata')
    if ret is not None:
        login, account, password = ret
        login = login.strip('"')
        password = password.strip('"')
        auth = (login, password)
    else:
        print('Error: could not load earthdata credentials')
        return()

    if url_base is None:
        url_base = '{}'.format(ac.config['tact_thredds_url_merra2'])
        print('Using base URL {}'.format(url_base))

    if obase is None:
        obase = os.path.abspath(ac.config['grid_dir']) + '/merra2/'

    ## parse date
    dt = dateutil.parser.parse(isotime)
    isodate = dt.isoformat()[0:10]
    c_time = dt.hour + dt.minute/60 + dt.second/3600

    ## set up url
    url = '{}/{}/{}/MERRA2_400.inst3_3d_asm_Np.{}.nc4'.format(url_base, dt.year, str(dt.month).zfill(2), dt.strftime('%Y%m%d'))
    print(url)

    ## load dimensions
    parstr = 'lat[0:1:360],lev[0:1:41],lon[0:1:575],time[0:1:7]'
    url_ascii = '{}.ascii?{}'.format(url, parstr)
    ret = requests.get(url_ascii, auth = auth)
    s = ret.content.decode('utf-8').split('\n')
    dims = {}
    for line in s:
        sp = [t.strip() for t in line.split(',')]
        if sp[0] not in ['lat','lev','lon','time']: continue
        dims[sp[0]] = np.asarray(sp[1:]).astype(float)
    dims['hour'] = dims['time']/(60)

    ## get levels
    levels = dims['lev']
    nlevels = len(levels)

    ## determine steps
    lon_step = dims['lon'][1]-dims['lon'][0]
    lat_step = dims['lat'][1]-dims['lat'][0]
    time_step = dims['hour'][1]-dims['hour'][0]

    ## determine bounding cells
    ## time
    bound_t0 = c_time - (c_time % time_step)
    bound_t1 = c_time - (c_time % time_step) + time_step
    time_cells = [bound_t0, bound_t1]

    ## compute bounding grid cells
    bound_s = limit[0] - (limit[0] % lat_step)
    bound_n = limit[2] - (limit[2] % lat_step) + lat_step
    bound_w = limit[1] - (limit[1] % lon_step)
    bound_e = limit[3] - (limit[3] % lon_step) + lon_step

    ## compute lats needed
    lat_range = bound_n - bound_s
    nlat = int(1+(lat_range / lat_step))
    lat_cells = [bound_s + i*lat_step for i in range(nlat)]

    ## compute lons needed
    lon_range = bound_e - bound_w
    nlon = int(1+(lon_range / lon_step))
    lon_cells = [bound_w + i*lon_step for i in range(nlon)]

    ## offset west longitudes
    if dims['lon'][0] >= 0:
        lon_cells = [c_lon if c_lon >= 0 else c_lon + 360 for c_lon in lon_cells]

    ## bounding cell indices
    xbounds = [np.argsort(np.abs(dims['lon']-min(lon_cells)))[0], np.argsort(np.abs(dims['lon']-max(lon_cells)))[0]]
    ybounds = [np.argsort(np.abs(dims['lat']-min(lat_cells)))[0], np.argsort(np.abs(dims['lat']-max(lat_cells)))[0]]
    tbounds = [np.argsort(np.abs(dims['hour']-min(time_cells)))[0], np.argsort(np.abs(dims['hour']-max(time_cells)))[0]]
    lbounds = [0, nlevels-1]


    ##[ time= 0 ..7] [ lev= 0 ..41] [ lat= 0 ..360] [ lon= 0 ..575]
    ## get parameter string
    parstr = 'RH[{}:1:{}][{}:1:{}][{}:1:{}][{}:1:{}]'.format(tbounds[0], tbounds[1], lbounds[0], lbounds[1],
                                                             ybounds[0], ybounds[1], xbounds[0], xbounds[1])
    parstr += ','
    parstr += 'T[{}:1:{}][{}:1:{}][{}:1:{}][{}:1:{}]'.format(tbounds[0], tbounds[1], lbounds[0], lbounds[1],
                                                             ybounds[0], ybounds[1], xbounds[0], xbounds[1])

    ## get data
    url_ascii = '{}.ascii?{}'.format(url, parstr)
    print(url_ascii)
    ret = requests.get(url_ascii, auth = auth)
    s = ret.content.decode('utf-8').split('\n')

    ## parse data
    pars = ['RH', 'T']
    data = {par: np.zeros((len(time_cells), nlevels, len(lat_cells), len(lon_cells))) + np.nan for par in pars}
    for l in s:
        if 'Dataset' in l: continue
        if len(l) == 0: continue
        sp = [t.strip() for t in l.split(',')]

        idx = [t.strip(']') for t in sp[0].split('[')]
        par = idx[0]
        idx = [int(i) for i in idx[1:]]
        for vi, v in enumerate(sp[1:]):
            v0 = float(v)
            if v0 < 1e15: data[par][idx[0], idx[1], idx[2], vi] = v0

    ## convert RH from fraction to percent
    data['RH'] *= 100

    ## check whether these files exist
    new_data = True if override else False
    if not new_data:
        for i,la in enumerate(lat_cells):
            for j,lo in enumerate(lon_cells):
                for k,ti in enumerate(time_cells):
                    odir = '{}/{}/{}/{}/{}'.format(obase, isodate, ti, la, lo)
                    if not os.path.exists(odir): os.makedirs(odir)
                    for par in pars:
                        ofile = '{}/{}.json'.format(odir, '_'.join([str(s) for s in [isodate, ti, la, lo, par]]))
                        if not os.path.exists(ofile):
                            new_data = True

    ## save profiles
    if verbosity > 0: print('Saving individual profiles')
    for i,la in enumerate(lat_cells):
            for j,lo in enumerate(lon_cells):
                for k,ti in enumerate(time_cells):
                    odir = '{}/{}/{}/{}/{}'.format(obase, isodate, ti, la, lo)
                    if not os.path.exists(odir): os.makedirs(odir)
                    for par in pars:
                        ofile = '{}/{}.json'.format(odir, '_'.join([str(s) for s in [isodate, ti, la, lo, par]]))

                        ## remove levels with NaNs
                        cur_data = [float(s) for s in list(data[par][k,:,i,j].data)]
                        cur_levels = [float(l) for l in levels]
                        cur_idx = [ind for ind, v in enumerate(cur_data) if np.isfinite(v)]
                        cur_data = [cur_data[ind] for ind in cur_idx]
                        cur_levels = [cur_levels[ind] for ind in cur_idx]

                        res = {'time':float(ti),'levels':cur_levels,'lat':la, 'lon':lo, 'data':cur_data}
                        if (not os.path.exists(ofile)) or (override):
                            with open(ofile, 'w') as f:
                                f.write(json.dumps(res))

    ## read & reformat profiles
    to_run = []
    for i,la in enumerate(lat_cells):
        for j,lo in enumerate(lon_cells):
            for k,ti in enumerate(time_cells):
                odir = '{}/{}/{}/{}/{}'.format(obase, isodate, ti, la, lo)
                if not os.path.exists(odir): os.makedirs(odir)

                ## write reformatted profile data
                tmp_profile = '{}/reformatted.profile'.format(odir)

                ## tmp to remove bad profiles
                if os.path.exists(tmp_profile):
                    os.remove(tmp_profile)

                if (not os.path.exists(tmp_profile)) or (override):
                    prof = {}
                    for par in ('RH', 'T'):
                        ofile = '{}/{}.json'.format(odir, '_'.join([str(s) for s in [isodate, ti, la, lo, par]]))
                        if os.path.exists(ofile):
                            prof[par] = json.load(open(ofile, 'r'))

                            ## invert profile
                            prof[par]['data'].reverse()
                            prof[par]['levels'].reverse()

                    ## we have the profiles
                    if len(prof)==2:
                        if verbosity > 1: print('Reformatting profiles {}'.format(tmp_profile))
                        with open(tmp_profile, 'w') as f:
                            f.write('{}\n'.format('# converted from MERRA2 profile'))
                            f.write('{}\n'.format('#   p(hPa)  T(K)  h2o(relative humidity)'))

                            for sj in range(len(prof['RH']['levels'])):
                                si = sj+0
                                f.write('{}\n'.format(' '.join([str(s) for s in (prof['RH']['levels'][si],
                                                                                 max(0,prof['T']['data'][si]),
                                                                                 max(0,prof['RH']['data'][si]))])))
                ## do simulation or read outputs
                if os.path.exists(tmp_profile):
                    to_run.append(tmp_profile)

    return((lat_cells, lon_cells, time_cells), to_run)
