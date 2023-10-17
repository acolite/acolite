# def tact_limit
## runs tact for a given limit file
## written by Quinten Vanhellemont, RBINS
## 2019-10-02
## modifications: 2019-12-17 renamed, added tact config, and removed dependencies
##                2021-02-27 (QV) integrated in acolite, added interpolation for target lat lon
##                2022-02-15 (QV) added L9/TIRS
##                2022-08-02 (QV) moved era5 profiles to separate function, added gdas1 profile option
##                2022-08-11 (QV) added ecostress and reptran

def tact_limit(isotime, limit=None,
                  lat = None, lon = None, wave_range=[7, 14],
                  source = 'era5', reptran = 'medium',
                  satsen = None, override = False, verbosity = 0, processes = 4):

    import netCDF4
    import numpy as np
    import scipy.interpolate
    import os, json, glob
    from functools import partial
    import multiprocessing
    import acolite as ac
    import datetime, dateutil.parser

    ## get verbosity from run settings
    verbosity = ac.settings['run']['verbosity']

    dt = dateutil.parser.parse(isotime)
    isodate = dt.isoformat()[0:10]
    c_time = dt.hour + dt.minute/60 + dt.second/3600

    source_default = 'era5'
    if source not in ['era5', 'gdas1', 'ncep.reanalysis2', 'merra2']: # 'ncep.reanalysis',
        print('Source {} not configured. Using {}.'.format(source, source_default))
        source = '{}'.format(source_default)

    if limit is None:
        if (lat is None) and (lon is None):
            return()
        limit = [np.nanmin(lat), np.nanmin(lon),
                 np.nanmax(lat), np.nanmax(lon)]

    ## read thermal band RSR
    rsr_data = {}
    for sen in ac.config['thermal_sensors']:
        if sen in ['L5_TM', 'L7_ETM']:
            rsr_file = "{}/RSR/{}_B6.txt".format(ac.config['data_dir'], sen)
        else:
            rsr_file = "{}/RSR/{}.txt".format(ac.config['data_dir'], sen)
        r_, b_ = ac.shared.rsr_read(rsr_file)
        rsr_data[sen]={'rsr':r_, 'bands':b_}

    ## directory to store retrieved profiles and simulation results
    obase = os.path.abspath(ac.config['grid_dir']) + '/{}/'.format(source)
    if source == 'era5':
        cells, to_run = ac.tact.tact_profiles_era5(isotime, limit, obase = obase, override = override, verbosity = verbosity)
    if source == 'gdas1':
        cells, to_run = ac.tact.tact_profiles_gdas1(isotime, limit, obase = obase, override = override, verbosity = verbosity)
    if source in ['merra2']:
        cells, to_run = ac.tact.tact_profiles_merra2(isotime, limit, obase = obase, override = override, verbosity = verbosity)
    if source in ['ncep.reanalysis', 'ncep.reanalysis2']:
        cells, to_run = ac.tact.tact_profiles_ncep(isotime, limit, obase = obase, override = override, verbosity = verbosity, source = source)

    ## space/time cells
    lat_cells, lon_cells, time_cells = cells

    if verbosity > 1: print('Running simulations for TACT using libRadtran at {}'.format(ac.config['libradtran_dir']))
    ## run stuff in multiprocessing
    with multiprocessing.Pool(processes=processes) as pool:
        results = pool.map(partial(ac.tact.tact_simulations, atmosphere=None, reptran = reptran, wave_range=wave_range,
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
