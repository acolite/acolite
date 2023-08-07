# def tact_profiles_ncep
## gets and reformats ncep reanalysis profiles for a given limit and dataset
##
## more info on this dataset at NOAA PSL
## https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html
## https://psl.noaa.gov/data/gridded/data.ncep.reanalysis2.html
##
## Don't use ncep.reanalysis as the profile levels are only between 300 and 1000 hPa
##
## written by Quinten Vanhellemont, RBINS
## 2022-08-03
## modifications: 2023-08-07 (QV) moved url_base to ACOLITE config file
##

def tact_profiles_ncep(isotime, limit, obase = None, override = False, verbosity = 5,
              source = 'ncep.reanalysis2', url_base = None, geo_step = 2.5, time_step = 6):

    import os
    import netCDF4
    import json
    import datetime, dateutil.parser
    import dateutil
    import numpy as np
    import acolite as ac

    source_default = 'ncep.reanalysis2'
    if source not in ['ncep.reanalysis2']: # 'ncep.reanalysis',
        print('Source {} not configured. Using {}.'.format(source, source_default))
        source = '{}'.format(source_default)

    if obase is None:
        obase = os.path.abspath(ac.config['grid_dir']) + '/{}/'.format(source)

    if url_base is None:
        url_base = '{}'.format(ac.config['tact_thredds_url_ncep'])
        print('Using base URL {}'.format(url_base))

    ## parse date
    dt = dateutil.parser.parse(isotime)
    isodate = dt.isoformat()[0:10]
    c_time = dt.hour + dt.minute/60 + dt.second/3600

    ## dataset start time
    dst = datetime.datetime(1800,1,1, tzinfo=dateutil.tz.tzutc())

    ## current time in hours since 1800-1-1
    td = dt - datetime.datetime(1800,1,1, tzinfo=dateutil.tz.tzutc())
    ct = td.days*24 + td.seconds/3600

    ## compute bounding grid cells
    bound_s = limit[0] - (limit[0] % geo_step)
    bound_n = limit[2] - (limit[2] % geo_step) + geo_step
    bound_w = limit[1] - (limit[1] % geo_step)
    bound_e = limit[3] - (limit[3] % geo_step) + geo_step

    ## compute lats needed
    lat_range = bound_n - bound_s
    nlat = int(1+(lat_range / geo_step))
    lat_cells = [bound_s + i*geo_step for i in range(nlat)]

    ## compute lons needed
    lon_range = bound_e - bound_w
    nlon = int(1+(lon_range / geo_step))
    lon_cells = [bound_w + i*geo_step for i in range(nlon)]

    ## offset west longitudes
    lon_cells = [c_lon if c_lon >= 0 else c_lon + 360 for c_lon in lon_cells]

    ## bounding time steps
    mt = np.mod(ct, time_step)
    time_cells = [ct - mt, ct + (time_step-mt)]

    for par in ['r', 't']:
        if par == 'r':
            dsi = 'rhum'
        elif par == 't':
            dsi = 'air'

        ## construct url
        url = '{}/{}/pressure/{}.{}.nc'.format(url_base, source, dsi, dt.year)

        ## check whether these files exist
        new_data = True if override else False
        if not new_data:
            for i,la in enumerate(lat_cells):
                for j,lo in enumerate(lon_cells):
                    for k,ti in enumerate(time_cells):
                        cur_time = dst + datetime.timedelta(seconds=ti*3600)
                        isodate = cur_time.isoformat()[0:10]
                        cur_hour = cur_time.hour

                        odir = '{}/{}/{}/{}/{}'.format(obase, isodate, cur_hour, la, lo)
                        if not os.path.exists(odir): os.makedirs(odir)
                        ofile = '{}/{}.json'.format(odir, '_'.join([str(s) for s in [isodate, cur_hour, la, lo, par]]))
                        if not os.path.exists(ofile):
                            new_data = True

        ## get new data
        if new_data:
            if verbosity > 0:
                print('Opening NetCDF {}'.format(url))

            ## open NetCDF
            ds = netCDF4.Dataset(url)
            datasets = ds.variables.keys()

            ## get dimensions
            times = ds.variables['time'][:]
            levels = ds.variables['level'][:]
            lats = ds.variables['lat'][:]
            lons = ds.variables['lon'][:]

            ## get geo and time bounds
            xbounds = [np.argsort(np.abs(lons-min(lon_cells)))[0], np.argsort(np.abs(lons-max(lon_cells)))[0]]
            ybounds = [np.argsort(np.abs(lats-min(lat_cells)))[0], np.argsort(np.abs(lats-max(lat_cells)))[0]]
            tidx = np.argsort(np.abs(times-ct))[0:2]

            ## get lat and lon steps (should be == to lat_cells and lon_cells)
            blon = lons[xbounds[0]:xbounds[1]].data
            blat = lats[ybounds[1]:ybounds[0]].data
            btime = times[tidx].data

            ## get bounding profiles for the four locations at times tidx
            prof = ds.variables[dsi][tidx, :, ybounds[1]:ybounds[0]+1, xbounds[0]:xbounds[1]+1]
            ds.close()
            if verbosity > 0: print('Closed NetCDF {}'.format(url))

            ## save profiles
            if verbosity > 0: print('Saving individual profiles')

            for i,la in enumerate(lat_cells):
                for j,lo in enumerate(lon_cells):
                    for k,ti in enumerate(time_cells):
                        cur_time = dst + datetime.timedelta(seconds=ti*3600)
                        isodate = cur_time.isoformat()[0:10]
                        cur_hour = cur_time.hour

                        odir = '{}/{}/{}/{}/{}'.format(obase, isodate, cur_hour, la, lo)
                        if not os.path.exists(odir): os.makedirs(odir)
                        ofile = '{}/{}.json'.format(odir, '_'.join([str(s) for s in [isodate, cur_hour, la, lo, par]]))

                        res = {'time':float(cur_hour),
                               'levels':[float(s) for s in list(levels.data)],
                               'lat':la, 'lon':lo,
                               'data':[float(s) for s in list(prof[k, :, len(lat_cells)-1-i, j].data)]}
                        if res['levels'][0] > res['levels'][-1]:
                            res['levels'].reverse()
                            res['data'].reverse()

                        if (not os.path.exists(ofile)) or (override):
                            with open(ofile, 'w') as f:
                                f.write(json.dumps(res))

    ## read & reformat profiles
    to_run = []
    for i,la in enumerate(lat_cells):
        for j,lo in enumerate(lon_cells):
            for k,ti in enumerate(time_cells):
                cur_time = dst + datetime.timedelta(seconds=ti*3600)
                isodate = cur_time.isoformat()[0:10]
                cur_hour = cur_time.hour

                odir = '{}/{}/{}/{}/{}'.format(obase, isodate, cur_hour, la, lo)
                if not os.path.exists(odir): os.makedirs(odir)

                ## write reformatted profile data
                tmp_profile = '{}/reformatted.profile'.format(odir)

                ## tmp to remove bad profiles
                if os.path.exists(tmp_profile):
                    os.remove(tmp_profile)

                if (not os.path.exists(tmp_profile)) or (override):
                    prof = {}
                    for par in ('r', 't'):
                        ofile = '{}/{}.json'.format(odir, '_'.join([str(s) for s in [isodate, cur_hour, la, lo, par]]))
                        if os.path.exists(ofile):
                            prof[par] = json.load(open(ofile, 'r'))

                    ## we have the profiles
                    if len(prof)==2:
                        if verbosity > 1: print('Reformatting profiles {}'.format(tmp_profile))
                        with open(tmp_profile, 'w') as f:
                            f.write('{}\n'.format('# converted from NCEP reanalysis profile'))
                            f.write('{}\n'.format('#   p(hPa)  T(K)  h2o(relative humidity)'))

                            for sj in range(len(prof['r']['levels'])):
                                si = sj+0
                                f.write('{}\n'.format(' '.join([str(s) for s in (prof['r']['levels'][si],
                                                                                 max(0,prof['t']['data'][si]),
                                                                                 max(0,prof['r']['data'][si]))])))

                ## do simulation or read outputs
                if os.path.exists(tmp_profile):
                    to_run.append(tmp_profile)

    time_cells_hour = []
    for k,ti in enumerate(time_cells):
        cur_time = dst + datetime.timedelta(seconds=ti*3600)
        isodate = cur_time.isoformat()[0:10]
        time_cells_hour.append(cur_time.hour)

    return((lat_cells, lon_cells, time_cells_hour), to_run)
