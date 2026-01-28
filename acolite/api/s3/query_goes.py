## def query_goes
## simple query and download function for GOES on AWS S3 storage
## written by Quinten Vanhellemont, RBINS
## 2025-12-17
## modifications: 2026-01-28 (QV) added sza threshold option

def query_goes(start_date, end_date = None, local_directory = None, time_range_sec = 600, download = False, override = False,
               satellite_index = 19, channels = ['M6C01', 'M6C02', 'M6C03', 'M6C04', 'M6C05', 'M6C06'],
               max_sun_zenith_angle = 70, station_lon = None, station_lat = None,
               url_base = 's3://noaa-goes', product_base = 'ABI-L1b-Rad', scan_factor = 'C', verbosity = 0, ):

    import os, dateutil.parser, datetime
    import s3fs
    import acolite as ac

    ## compute sun zenith angles for location
    compute_sun = (max_sun_zenith_angle is not None) &\
                    (station_lon is not None) & (station_lat is not None)

    ## check scan factor
    scan_factors = ['C', 'F', 'M1', 'M2']
    if scan_factor not in scan_factors:
        print('Scan factor {} not configured. Use one of {}')
        return
    product = '{}{}'.format(product_base, scan_factor)

    ## parse dates
    dt_start = dateutil.parser.parse(start_date)
    if end_date is None:
        dt_end = dt_start + datetime.timedelta(seconds = time_range_sec)
    else:
        dt_end = dateutil.parser.parse(end_date)

    dates = []
    dt_cur = datetime.datetime(dt_start.year, dt_start.month, dt_start.day)
    dt_max = datetime.datetime(dt_end.year, dt_end.month, dt_end.day, 23, 59, 59)
    while dt_cur < dt_max:
        dates.append({'datetime': dt_cur, 'isodate': dt_cur.isoformat()[0:10], 'doy': dt_cur.strftime('%j').zfill(3), 'year': str(dt_cur.year)})
        dt_cur += datetime.timedelta(days = 1)
    if len(dates) == 0: return

    if local_directory is None:
        local_directory = os.getcwd()
        if verbosity > 0: print('Warning: Downloading to current directory: {}'.format(local_directory))

    ## set up anonymous S3 file system
    fs = s3fs.S3FileSystem(anon = True)

    remote_files = []
    local_files = []

    ## run through dates
    for di, dt in enumerate(dates):
        files = fs.glob('{}{}/{}/{}/{}/*/*'.format(url_base, satellite_index, product, dt['year'], dt['doy']))
        if verbosity > 1: print('Found {} total files for date {}'.format(len(files), dt['isodate']))
        n0 = len(remote_files)

        ## run through files
        for file in files:
            ## check local file
            local_file = '{}/{}'.format(local_directory, file)
            if os.path.exists(local_file):
                if (override): os.remove(local_file)
            if not os.path.exists(local_file):
                ## parse file name
                bn = os.path.basename(file)
                bn_ = bn[bn.find(product)+len(product)+1:-3]
                sp = bn_.split('_')

                ## skip if not in listed bands
                if sp[0] not in channels: continue

                ## parse image date date
                sdate = sp[2][1:]
                h, m, s = sdate[7:9], sdate[9:11], sdate[11:13]
                image_isodate = dt['isodate']+'T{}:{}:{}'.format(h,m,s)
                image_dtime = dateutil.parser.parse(image_isodate)

                ## compute sun position for this image and given location
                if compute_sun:
                    sun_angle = ac.shared.sun_position(image_dtime, float(station_lon), float(station_lat))['zenith']
                    if sun_angle > float(max_sun_zenith_angle): continue ## sun too low

                if (image_dtime < dt_start): continue ## continue if before start date
                if (image_dtime > dt_end): continue ## continue if after end date

                ## download file
                if download:
                    if verbosity > 2: print('Downloading {}'.format(file))
                    fs.download(file,local_file)
                    if verbosity > 2: print('Downloaded to {}'.format(local_file))

            remote_files.append(file)
            if os.path.exists(local_file):
                local_files.append(local_file)

        if verbosity > 1: print('Selected {} files from date {}'.format(len(remote_files) - n0, dt['isodate']))

    if download:
        return(local_files)
    else:
        return(remote_files)
