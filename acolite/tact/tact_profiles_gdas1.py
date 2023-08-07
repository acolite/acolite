# def tact_profiles_gdas1
## gets and reformats gdas1 profiles for a given limit and dataset
##
## more info on this dataset at RDA UCAR
## https://rda.ucar.edu/datasets/ds083.3/
## ds083.3 | DOI: 10.5065/D65Q4T4Z
##
## written by Quinten Vanhellemont, RBINS
## 2022-08-02
## modifications: 2023-08-07 (QV) moved url_base to ACOLITE config file

def tact_profiles_gdas1(isotime, limit, obase = None, override = False, verbosity = 5,
               url_base = None, geo_step = 0.25):

    import os
    import netCDF4
    import json
    import dateutil.parser
    import numpy as np
    import acolite as ac
    import datetime

    if obase is None:
        os.path.abspath(ac.config['grid_dir']) + '/gdas1/'

    if url_base is None:
        url_base = '{}'.format(ac.config['tact_thredds_url_gdas1'])
        print('Using base URL {}'.format(url_base))

    ## parse date
    dt = dateutil.parser.parse(isotime)
    isodate = dt.isoformat()[0:10]
    c_time = dt.hour + dt.minute/60 + dt.second/3600

    ## compute times
    time_step = 6 * 3600
    times = np.array([0, 6, 12, 18])
    for ti, time in enumerate(times):
        if c_time < time: continue
        idx = ti

    ## bounding times
    d0 = datetime.datetime(dt.year, dt.month, dt.day, hour=times[idx])
    d1 = d0 + datetime.timedelta(seconds=time_step)

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

    ## time cells
    #time_cells = [times[idx], 24 if idx == 3 else times[idx+1]]
    time_cells = [d0.hour, d1.hour]

    ## check whether these files exist
    new_data = True if override else False
    if not new_data:
        for i,la in enumerate(lat_cells):
            for j,lo in enumerate(lon_cells):
                for k,ti in enumerate(time_cells):
                    odir = '{}/{}/{}/{}/{}'.format(obase, isodate, ti, la, lo)
                    if not os.path.exists(odir): os.makedirs(odir)
                    for par in ['r', 't']:
                        ofile = '{}/{}.json'.format(odir, '_'.join([str(s) for s in [isodate, ti, la, lo, par]]))
                        if not os.path.exists(ofile):
                            new_data = True

    if new_data:
        ## run through time steps
        for d in [d0, d1]:
                btime = d.hour
                date = d.strftime('%Y%m%d%H')
                url = '{}/{}/{}/gdas1.fnl0p25.{}.f00.grib2'\
                            .format(url_base, date[0:4], date[0:6], date)

                if verbosity > 0:
                    print('Opening NetCDF {}'.format(url))

                prof = {}

                ## open NetCDF
                ds = netCDF4.Dataset(url)
                datasets = ds.variables.keys()

                ## run through datasets
                for dsi in ['Relative_humidity_isobaric', 'Temperature_isobaric']:
                    k = dsi[0].lower()
                    ## get dimensions
                    #levels = ds['isobaric'][:]
                    coords = ds[dsi].getncattr('coordinates').split()
                    levels = ds[coords[-3]][:]
                    print(coords[-3], len(levels))
                    lats = ds['lat'][:]
                    lons = ds['lon'][:]
                    times = [d.hour]

                    ## get geo and time bounds
                    xbounds = [np.argsort(np.abs(lons-min(lon_cells)))[0], np.argsort(np.abs(lons-max(lon_cells)))[0]]
                    ybounds = [np.argsort(np.abs(lats-min(lat_cells)))[0], np.argsort(np.abs(lats-max(lat_cells)))[0]]
                    tidx = np.argsort(np.abs(time-c_time))[0:2]

                    ## get lat and lon steps (should be == to lat_cells and lon_cells)
                    blon = lons[xbounds[0]:xbounds[1]].data
                    blat = lats[ybounds[1]:ybounds[0]].data

                    ## get bounding profiles for the four locations at times tidx
                    prof[k] = {}
                    prof[k]['levels'] = levels
                    prof[k]['data'] = ds.variables[dsi][0, :, ybounds[1]:ybounds[0]+1, xbounds[0]:xbounds[1]+1]
                ## close NetCDF
                ds.close()
                if verbosity > 0: print('Closed NetCDF {}'.format(url))

                ## save profiles
                if verbosity > 0: print('Saving individual profiles')
                ti = '{}'.format(d.hour)
                for par in prof:
                    for i,la in enumerate(lat_cells):
                        for j,lo in enumerate(lon_cells):
                            odir = '{}/{}/{}/{}/{}'.format(obase, isodate, ti, la, lo)
                            if not os.path.exists(odir): os.makedirs(odir)
                            ofile = '{}/{}.json'.format(odir, '_'.join([str(s) for s in [isodate, ti, la, lo, par]]))

                            res = {'time':float(ti),
                                   'levels':[float(s)/100 for s in prof[par]['levels']],
                                   'lat':la, 'lon':lo,
                                   #'data':[float(s) for s in list(prof[par]['data'][:, i, j].data)]
                                   'data':[float(s) for s in list(prof[par]['data'][:, len(lat_cells)-1-i, j].data)]
                                  }

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
                    for par in ('r', 't'):
                        ofile = '{}/{}.json'.format(odir, '_'.join([str(s) for s in [isodate, ti, la, lo, par]]))
                        if os.path.exists(ofile):
                            prof[par] = json.load(open(ofile, 'r'))

                    ## we have the profiles
                    if len(prof)==2:
                        if verbosity > 1: print('Reformatting profiles {}'.format(tmp_profile))
                        with open(tmp_profile, 'w') as f:
                            f.write('{}\n'.format('# converted from GDAS1 profile'))
                            f.write('{}\n'.format('#   p(hPa)  T(K)  h2o(relative humidity)'))

                            for sj in range(len(prof['r']['levels'])):
                                si = sj+0
                                f.write('{}\n'.format(' '.join([str(s) for s in (prof['r']['levels'][si],
                                                                                 max(0,prof['t']['data'][si]),
                                                                                 max(0,prof['r']['data'][si]))])))

                ## do simulation or read outputs
                if os.path.exists(tmp_profile):
                    to_run.append(tmp_profile)


    return((lat_cells, lon_cells, time_cells), to_run)
