## def tact.profiles.era5_ecmwf
## gets and reformats era5 profiles for a given limit and dataset
## needs conda install -c conda-forge ecmwf-datastores-client
##
## here the data are retrieved using ecmwf.datastores
## https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels?tab=overview
##
## written by Quinten Vanhellemont, RBINS
## 2026-05-26
## modifications: 2026-05-26 (QV) added to tact.profiles

def era5_ecmwf(isotime, limit, obase = None, override = False, verbosity = 5, geo_step = 0.25, delete = True):
    import os, time
    import netCDF4
    import json
    import datetime, dateutil.parser
    import numpy as np
    import acolite as ac

    if obase is None: obase = os.path.abspath(ac.config['grid_dir']) + '/era5_ecmwf'

    ## parse date
    dt = dateutil.parser.parse(isotime)
    isodate = dt.isoformat()[0:10]
    c_time = dt.hour + dt.minute/60 + dt.second/3600

    ## compute bounding grid cells
    bound_s = limit[0] - (limit[0] % geo_step)
    bound_n = limit[2] - (limit[2] % geo_step) + geo_step
    bound_w = limit[1] - (limit[1] % geo_step)
    bound_e = limit[3] - (limit[3] % geo_step) + geo_step

    ## compute lats needed
    lat_range = bound_n - bound_s
    nlat = int(1+(lat_range / geo_step))
    lat_cells = [bound_s + i*geo_step for i in range(nlat)]
    lat_cells.reverse() ## sort from North to South

    ## compute lons needed
    lon_range = bound_e - bound_w
    nlon = int(1+(lon_range / geo_step))
    lon_cells = [bound_w + i*geo_step for i in range(nlon)]

    ## offset west longitudes
    lon_cells = [c_lon if c_lon >= 0 else c_lon + 360 for c_lon in lon_cells]

    ## time cells
    time_cells = [int(np.floor(c_time)), int(np.floor(c_time))+1]

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

    ## set up query
    dataset = "reanalysis-era5-pressure-levels"
    request = {
        "product_type": ["reanalysis"],
        "variable": ["relative_humidity", "temperature"],
        "year": [str(dt.year)], "month": [str(dt.month).zfill(2)], "day": [str(dt.day).zfill(2)],
        "time": [str(time_cells[0]).zfill(2) + ":00", str(time_cells[1]).zfill(2) + ":00"],
        "pressure_level": ["1", "2", "3", "5", "7", "10","20", "30", "50", "70", "100", "125", "150", "175",
                           "200", "225", "250", "300", "350", "400", "450", "500", "550", "600", "650", "700",
                           "750", "775", "800", "825", "850", "875", "900", "925", "950", "975", "1000"],
        "data_format": "netcdf", "download_format": "unarchived",
        ## data are always sorted from N->S W->E
        "area": [str(v) for v in [lat_cells[0], lon_cells[0], lat_cells[-1], lon_cells[-1]]],
        }

    ## temporary name
    base_name = 'ERA5_profiles_rt_{}_{}_{}_{}_{}.nc'.format(request['year'][0], request['month'][0], request['day'][0],
                                                            '_'.join(request['time']), '_'.join(request['area']))
    local_netcdf = ac.config['scratch_dir'] + '/' + base_name

    ## download data
    if new_data:
        if not os.path.exists(local_netcdf) | override:
            from ecmwf.datastores import Client
            client = Client()

            ## verify connection
            try:
                client.check_authentication()
            except BaseException as err:
                print("Download error {}, {}".format(err, type(err)))
                print('Please configure url and key in .ecmwfdatastoresrc')
                return

            ## retrieve data
            client.retrieve(dataset, request, target = local_netcdf)

        ## load data
        if not os.path.exists(local_netcdf):
            print('Could not get required data from ECMWF.')
            return

        ## get grid info
        lat = ac.shared.nc_data(local_netcdf, 'latitude')
        lon = ac.shared.nc_data(local_netcdf, 'longitude')
        times = ac.shared.nc_data(local_netcdf, 'valid_time')
        pressure_level = ac.shared.nc_data(local_netcdf, 'pressure_level')
        ## flip levels to store data from low->high pressure
        flip = True if pressure_level.data[0] > pressure_level.data[-1] else False
        if flip: levels = np.flip(pressure_level)

        ## run through parameters relative humidity (r) and temperature (t)
        for par in ['r', 't']:
            data = ac.shared.nc_data(local_netcdf, par)

            for k, ti in enumerate(time_cells):
                dt_ = datetime.datetime.fromtimestamp(times[k], tz=datetime.UTC)
                if (dt_.hour != ti): break

                for i,la in enumerate(lat_cells):
                    if (lat[i] != la): break

                    for j,lo in enumerate(lon_cells):
                        if (lon[j] != lo): break

                        odir = '{}/{}/{}/{}/{}'.format(obase, isodate, ti, la, lo)
                        if not os.path.exists(odir): os.makedirs(odir)
                        ofile = '{}/{}.json'.format(odir, '_'.join([str(s) for s in [isodate, ti, la, lo, par]]))

                        ## write result
                        if not os.path.exists(ofile) or (override):
                            if flip:
                                data_list = list(np.flip(data.data[k, :, i, j]))
                            else:
                                data_list = list(data.data[k, :, i, j])
                            res = {'time': float(ti),'levels':[float(l) for l in levels.data],
                                   'lat':float(la), 'lon':float(lo),
                                    'data':[float(s) for s in data_list]}
                            with open(ofile, 'w', encoding = 'utf-8') as f:
                                f.write(json.dumps(res))
                                
        ## delete downloaded netcdf
        if delete:
            os.remove(local_netcdf)
            print('Deleted {}'.format(local_netcdf))

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
                if os.path.exists(tmp_profile): os.remove(tmp_profile)

                if (not os.path.exists(tmp_profile)) or (override):
                    prof = {}
                    for par in ('r', 't'):
                        ofile = '{}/{}.json'.format(odir, '_'.join([str(s) for s in [isodate, ti, la, lo, par]]))
                        if os.path.exists(ofile):
                            prof[par] = json.load(open(ofile, 'r', encoding = 'utf-8'))

                    ## we have the profiles
                    if len(prof)==2:
                        if verbosity > 1: print('Reformatting profiles {}'.format(tmp_profile))
                        with open(tmp_profile, 'w', encoding = 'utf-8') as f:
                            f.write('{}\n'.format('# converted from ERA5 ECMWF profile'))
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
