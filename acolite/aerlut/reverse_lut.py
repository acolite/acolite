## read reverse luts and set up rgi
## luts are created if they do not exist
## reverse luts go from rpath -> aot
## QV 2021-02-03
## last updates: 2021-05-31 (QV) added remote lut retrieval
##               2021-10-24 (QV) added pressures and get_remote as keyword to other functions
##               2021-10-25 (QV) test if the wind dimension is != 1 or missing
##                2023-08-03 (QV) get lut url from ac.config

def reverse_lut(sensor, lutdw=None, par = 'romix',
                       pct = (1,60), nbins = 20, override = False,
                       pressures = [500, 750, 1013, 1100],
                       base_luts = ['ACOLITE-LUT-202110-MOD1', 'ACOLITE-LUT-202110-MOD2'],
                       rsky_lut = 'ACOLITE-RSKY-202102-82W',
                       get_remote = True, remote_base = None):
    import acolite as ac
    import numpy as np
    from netCDF4 import Dataset
    import scipy.interpolate
    import time, os

    ## use URL from main config
    if remote_base is None: remote_base = '{}'.format(ac.config['lut_url'])

    if lutdw is None:
        rsrf = ac.config['data_dir']+'/RSR/{}.txt'.format(sensor)
        rsr, rsr_bands = ac.shared.rsr_read(rsrf)
        bands = [b for b in rsr_bands]
    else:
        lut = list(lutdw.keys())[0]
        bands = list(lutdw[lut]['rgi'].keys())

    revl = {}
    for lut in base_luts:
        lutdir = '{}/{}-Reverse/{}'.format(ac.config['lut_dir'], '-'.join(lut.split('-')[0:3]), sensor)
        if not os.path.exists(lutdir): os.makedirs(lutdir)

        rgi = {}
        for b in bands:
            slut = '{}-reverse-{}-{}-{}'.format(lut, sensor, par, b)
            lutnc = '{}/{}.nc'.format(lutdir, slut)

            if (not os.path.exists(lutnc)) or (override):
                if os.path.exists(lutnc): os.remove(lutnc)

                ## try downloading LUT from GitHub
                if (get_remote):
                    remote_lut = '{}/{}-Reverse/{}/{}.nc'.format(remote_base, '-'.join(lut.split('-')[0:3]), sensor, slut)
                    try:
                        print('Getting remote LUT {}'.format(remote_lut))
                        ac.shared.download_file(remote_lut, lutnc)
                        print('Testing LUT {}'.format(lutnc))
                        lutb, meta = ac.shared.lutnc_import(lutnc) # test LUT
                    except:
                        print('Could not download remote lut {} to {}'.format(remote_lut, lutnc))
                        if os.path.exists(lutnc): os.remove(lutnc)

                ## generate LUT if download did not work
                if (not os.path.exists(lutnc)):
                    print('Creating reverse LUTs for {}'.format(sensor))
                    if lutdw is None:
                        print('Importing source LUTs')
                        lutdw = ac.aerlut.import_luts(sensor = sensor, base_luts = base_luts,
                                                        lut_par = [par], par = par, return_lut_array = True,
                                                        pressures = pressures, get_remote = get_remote,
                                                        add_rsky = (par == 'romix+rsky_t') | (par == 'romix+ffss_toa'), rsky_lut = rsky_lut)
                    pid = lutdw[lut]['ipd'][par]
                    if len(lutdw[lut]['dim']) == 7:
                        wind_dim = True
                        pressures, pids, raas, vzas, szas, winds, aots = lutdw[lut]['dim']
                    else:
                        pressures, pids, raas, vzas, szas, aots = lutdw[lut]['dim']
                        wind_dim = False
                        winds = np.atleast_1d(2)


                    print('Starting {}'.format(slut))
                    t0 = time.time()
                    tmp = lutdw[lut]['lut'][b][:,pid,:,:,:,:,:].flatten()
                    tmp = np.log(tmp)
                    prc = np.nanpercentile(tmp, pct)
                    h = np.histogram(tmp, bins=nbins, range=prc)
                    rpath_bins = np.exp(h[1])

                    ## set up dimensions for lut
                    lut_dimensions = ('pressure','raa','vza','sza','wind','rho')
                    dim = [pressures, raas, vzas, szas, winds, rpath_bins]
                    dims = [len(d) for d in dim]
                    luta = np.zeros(dims) + np.nan
                    ii = 0
                    ni = np.product(dims[:-1])
                    for pi, pressure in enumerate(pressures):
                        for ri, raa in enumerate(raas):
                            for vi, vza in enumerate(vzas):
                                for si, sza in enumerate(szas):
                                    for wi, wind in enumerate(winds):
                                        if wind_dim:
                                            ret = lutdw[lut]['rgi'][b]((pressure, pid, raa, vza, sza, wind, aots))
                                        else:
                                            ret = lutdw[lut]['rgi'][b]((pressure, pid, raa, vza, sza, aots))
                                        luta[pi, ri, vi, si, wi, :] = np.interp(rpath_bins, ret, aots)
                                        ii+=1
                            print('{} {:.1f}%'.format(b, (ii/ni)*100), end='\r')
                    print('\nResampling {} took {:.1f}s'.format(slut, time.time()-t0))

                    ## write this sensor band lut
                    if os.path.exists(lutnc): os.remove(lutnc)
                    nc = Dataset(lutnc, 'w')
                    ## set attributes
                    setattr(nc, 'base', slut)
                    setattr(nc, 'aermod', lut[-1])
                    setattr(nc, 'aots', aots)
                    setattr(nc, 'lut_dimensions', lut_dimensions)
                    for di, dn in enumerate(lut_dimensions):
                        ## set attribute
                        setattr(nc, dn, dim[di])
                        ## create dimensions
                        nc.createDimension(dn, len(dim[di]))
                    ## write lut
                    var = nc.createVariable('lut',np.float32,lut_dimensions)
                    var[:] = luta.astype(np.float32)
                    nc.close()

            ## read LUT and make rgi
            if os.path.exists(lutnc):
                nc = Dataset(lutnc)
                meta = {}
                for attr in nc.ncattrs():
                    attdata = getattr(nc,attr)
                    if isinstance(attdata,str): attdata = attdata.split(',')
                    meta[attr]=attdata
                lutb = nc.variables['lut'][:]
                nc.close()

                try:
                    minaot = np.nanmin(meta['aots'])
                    maxaot = np.nanmax(meta['aots'])
                except:
                    minaot = 0.001
                    maxaot = 5

                ## band specific interpolator
                if len(np.atleast_1d(meta['wind'])) == 1:
                    rgi[b] = scipy.interpolate.RegularGridInterpolator([meta[k] for k in meta['lut_dimensions'] if k not in ['wind']],
                                                                     lutb[:,:,:,:,0,:],bounds_error=False, fill_value=None)
                else:
                    rgi[b] = scipy.interpolate.RegularGridInterpolator([meta[k] for k in meta['lut_dimensions']],
                                                                     lutb,bounds_error=False, fill_value=None)
        revl[lut]={'rgi':rgi, 'minaot':minaot, 'maxaot':maxaot,
                   'model':int(lut[-1]), 'meta':meta}
    return(revl)
