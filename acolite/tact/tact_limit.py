# def tact_limit
## runs tact for a given limit file
## written by Quinten Vanhellemont, RBINS
## 2019-10-02
## modifications: 2019-12-17 renamed, added tact config, and removed dependencies
##                2021-02-27 (QV) integrated in acolite, added interpolation for target lat lon
##                2022-02-15 (QV) added L9/TIRS

def tact_limit(isotime, limit=None,
                  lat = None, lon = None,
                  url_base = 'https://rda.ucar.edu/thredds/dodsC/files/g/ds633.0/e5.oper.an.pl',
                  geo_step = 0.25,
                  satsens=['L9_TIRS', 'L8_TIRS', 'L5_TM', 'L7_ETM'],
                  satsen = 'L8_TIRS', override = False, verbosity = 0, processes = 4):

    import netCDF4
    import numpy as np
    import scipy.interpolate
    import os, json, glob
    from functools import partial
    import multiprocessing
    import acolite as ac
    import datetime, dateutil.parser

    dt = dateutil.parser.parse(isotime)
    isodate = dt.isoformat()[0:10]
    c_time = dt.hour + dt.minute/60 + dt.second/3600

    if limit is None:
        if (lat is None) and (lon is None):
            return()
        limit = [np.nanmin(lat), np.nanmin(lon),
                 np.nanmax(lat), np.nanmax(lon)]

    ## read thermal band RSR
    rsr_data = {}
    for sen in satsens:
        if 'TIRS' in sen:
            rsr_file = "{}/RSR/{}.txt".format(ac.config['data_dir'], sen)
        else:
            rsr_file = "{}/RSR/{}_B6.txt".format(ac.config['data_dir'], sen)
        r_, b_ = ac.shared.rsr_read(rsr_file)
        rsr_data[sen]={'rsr':r_, 'bands':b_}

    ## directory to store retrieved profiles and simulation results
    obase = os.path.abspath(ac.config['grid_dir'])

    ### read the data from netcdf
    for par in ('r', 't'):
        date = isodate.replace('-','')[0:8]
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
            ## open NetCDF
            ds = netCDF4.Dataset(url)
            datasets = ds.variables.keys()

            ## get dimensions
            utc_date = ds.variables['utc_date'][:]
            date = (utc_date/100).astype(int)
            time = (utc_date-date*100).astype(int)
            levels = ds.variables['level'][:]
            lats = ds.variables['latitude'][:]
            lons = ds.variables['longitude'][:]
            times = ds.variables['time'][:]

            ## get geo and time bounds
            xbounds = [np.argsort(np.abs(lons-min(lon_cells)))[0], np.argsort(np.abs(lons-max(lon_cells)))[0]]
            ybounds = [np.argsort(np.abs(lats-min(lat_cells)))[0], np.argsort(np.abs(lats-max(lat_cells)))[0]]
            tidx = np.argsort(np.abs(time-c_time))[0:2]

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

                        odir = '{}/{}/{}/{}/{}'.format(obase, isodate, ti, la, lo)
                        if not os.path.exists(odir): os.makedirs(odir)
                        ofile = '{}/{}.json'.format(odir, '_'.join([str(s) for s in [isodate, ti, la, lo, par]]))

                        res = {'time':float(ti),'levels':list(levels),'lat':la, 'lon':lo,
                               'data':[float(s) for s in list(prof[k, :, i, j].data)]}

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

    ## run stuff in multiprocessing
    with multiprocessing.Pool(processes=processes) as pool:
        results = pool.map(partial(ac.tact.tact_simulations, atmosphere=None,
                                                pdate='', rsr_data=rsr_data, obase=None), to_run)


    ## read simulation outputs
    sims = None
    for i,la in enumerate(lat_cells):
        for j,lo in enumerate(lon_cells):
            for k,ti in enumerate(time_cells):
                odir = '{}/{}/{}/{}/{}'.format(obase, isodate, ti, la, lo)
                if not os.path.exists(odir): os.makedirs(odir)

                ## write reformatted profile data
                tmp_profile = '{}/reformatted.profile'.format(odir)
                if verbosity > 1: print('Importing simulation {}'.format(tmp_profile))

                ## import sim result
                sdir = '{}/reformatted'.format(odir)
                sf = glob.glob('{}/*{}.txt'.format(sdir, satsen))
                sim = {}
                if len(sf) >= 1:
                    sf = sf[0]
                    dt = os.path.basename(sf).split('_')[2]
                    sd = ac.tact.read_sim(sf)

                    for ib, b in enumerate(rsr_data[satsen]['bands']):
                        sim['tau{}'.format(b)] = sd['tau'][ib]
                        sim['Lu{}'.format(b)] = sd['Lu'][ib]
                        sim['Ld{}'.format(b)] = sd['Ld'][ib]

                    if sims is None:
                        sims = {s: np.zeros((len(lat_cells),len(lon_cells), len(time_cells))) for s in sim}
                    for s in sim:
                        sims[s][i,j,k] = sim[s]


    ## time interpolation
    ## calculate the time weighting
    w = np.interp(c_time, time_cells, (0,1))
    simst = {s: sims[s][:,:, 1] * w + sims[s][:,:,0] * (1.-w) for s in sims}

    ## return lons between -180 and 180
    for il, lonv in enumerate(lon_cells):
        if lonv > 180: lon_cells[il]-=360

    ## interpolate to given lat lon
    if (lat is not None) & (lon is not None):
        ## interpolate simulation data
        x,y = np.meshgrid(lon_cells,lat_cells)
        xr = list(x.ravel())
        yr = list(y.ravel())
        thd = {k: scipy.interpolate.griddata((xr,yr), list(simst[k].ravel()),
                                             (lon, lat), method='linear') for k in simst.keys()}
        return(thd, simst, lon_cells, lat_cells)
    else:
        return(simst, lon_cells, lat_cells)
