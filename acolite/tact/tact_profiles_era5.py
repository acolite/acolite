# def tact_profiles_era5
## gets and reformats era5 profiles for a given limit and dataset
##
## more info on this dataset at RDA UCAR
## https://rda.ucar.edu/datasets/ds633.0/
## ds633.0 | DOI: 10.5065/BH6N-5N20
##
## written by Quinten Vanhellemont, RBINS
## 2019-10-02
## modifications: 2022-08-02 (QV) split off from tact_limit.py
##                                flipped lats subset
##                2023-01-29 (QV) added iteration in case NetCDF I/O error occurs
##                2023-01-30 (QV) added grib option
##                2023-08-07 (QV) moved url_base to ACOLITE config file

def tact_profiles_era5(isotime, limit, obase = None, override = False, verbosity = 5, grib = False,
              url_base = None, geo_step = 0.25):

    import os, time
    import netCDF4
    import json
    import dateutil.parser
    import numpy as np
    import acolite as ac

    if obase is None:
        obase = os.path.abspath(ac.config['grid_dir']) + '/era5/'

    if url_base is None:
        url_base = '{}'.format(ac.config['tact_thredds_url_era5'])
        print('Using base URL {}'.format(url_base))

    ## parse date
    dt = dateutil.parser.parse(isotime)
    isodate = dt.isoformat()[0:10]
    c_time = dt.hour + dt.minute/60 + dt.second/3600

    ### read the data from netcdf
    for par in ('r', 't'):
        date = isodate.replace('-','')[0:8]

        if grib:
            latpar = 'lat'
            lonpar = 'lon'
            levpar = 'isobaric'
            datpar = 'time'
            if par == 'r':
                dsi = 'Relative_humidity_isobaric'
                url = '{}/{}/e5.oper.an.pl.128_157_r.ll025sc.{}00_{}23.grb'.format(url_base, date[0:6], date, date)
            elif par == 't':
                dsi = 'Temperature_isobaric'
                url = '{}/{}/e5.oper.an.pl.128_130_t.ll025sc.{}00_{}23.grb'.format(url_base, date[0:6], date, date)
        else:
            latpar = 'latitude'
            lonpar = 'longitude'
            levpar = 'level'
            datpar = 'utc_date'
            if par == 'r':
                dsi = 'R'
                url = '{}/{}/e5.oper.an.pl.128_157_r.ll025sc.{}00_{}23.nc'.format(url_base, date[0:6], date, date)
            elif par == 't':
                dsi = 'T'
                url = '{}/{}/e5.oper.an.pl.128_130_t.ll025sc.{}00_{}23.nc'.format(url_base, date[0:6], date, date)

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
        time_cells = [int(np.floor(c_time)), int(np.ceil(c_time))]

        ## check whether these files exist
        new_data = True if override else False
        if not new_data:
            for i,la in enumerate(lat_cells):
                for j,lo in enumerate(lon_cells):
                    for k,ti in enumerate(time_cells):

                        odir = '{}/{}/{}/{}/{}'.format(obase, isodate, ti, la, lo)
                        if not os.path.exists(odir): os.makedirs(odir)
                        ofile = '{}/{}.json'.format(odir, '_'.join([str(s) for s in [isodate, ti, la, lo, par]]))
                        if not os.path.exists(ofile):
                            new_data = True

        ## get new data
        if new_data:
            if verbosity > 0:
                print('Opening NetCDF {}'.format(url))

            it = 0
            success = False
            while (it < 2) & (success == False):
                try:
                    ## open NetCDF
                    ds = netCDF4.Dataset(url)
                    datasets = ds.variables.keys()

                    ## get dimensions
                    utc_date = ds.variables[datpar][:]
                    if grib: utc_date = np.asarray([np.int32('{}{}'.format(date, '{:.0f}'.format(t).zfill(2))) for t in utc_date])
                    date = (utc_date/100).astype(int)
                    ttime = (utc_date-date*100).astype(int)
                    levels = ds.variables[levpar][:]
                    lats = ds.variables[latpar][:]
                    lons = ds.variables[lonpar][:]
                    #times = ds.variables['time'][:]

                    ## get geo and time bounds
                    xbounds = [np.argsort(np.abs(lons-min(lon_cells)))[0], np.argsort(np.abs(lons-max(lon_cells)))[0]]
                    ybounds = [np.argsort(np.abs(lats-min(lat_cells)))[0], np.argsort(np.abs(lats-max(lat_cells)))[0]]
                    tidx = np.argsort(np.abs(ttime-c_time))[0:2]

                    ## get lat and lon steps (should be == to lat_cells and lon_cells)
                    blon = lons[xbounds[0]:xbounds[1]].data
                    blat = lats[ybounds[1]:ybounds[0]].data
                    #btime = times[tidx].data

                    ## get bounding profiles for the four locations at times tidx
                    prof = ds.variables[dsi][tidx, :, ybounds[1]:ybounds[0]+1, xbounds[0]:xbounds[1]+1]
                    ds.close()
                    if verbosity > 0: print('Closed NetCDF {}'.format(url))
                    success = True

                except BaseException as err:
                    print('Could not open {}'.format(url))
                    print("Error {}, {}".format(err, type(err)))
                    if it < 2: time.sleep(10)
                    pass
                it +=1
            if not success:
                print('Failed to open {}'.format(url))
                return()

            ## save profiles
            if verbosity > 0: print('Saving individual profiles')

            for i,la in enumerate(lat_cells):
                for j,lo in enumerate(lon_cells):
                    for k,ti in enumerate(time_cells):

                        odir = '{}/{}/{}/{}/{}'.format(obase, isodate, ti, la, lo)
                        if not os.path.exists(odir): os.makedirs(odir)
                        ofile = '{}/{}.json'.format(odir, '_'.join([str(s) for s in [isodate, ti, la, lo, par]]))

                        res = {'time':float(ti),'levels':[float(l) for l in levels],'lat':la, 'lon':lo,
                               #'data':[float(s) for s in list(prof[k, :, i, j].data)]
                               'data':[float(s) for s in list(prof[k, :, len(lat_cells)-1-i, j].data)]}
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
                            f.write('{}\n'.format('# converted from ERA5 profile'))
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
