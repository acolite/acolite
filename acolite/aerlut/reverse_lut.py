## read reverse luts and set up rgi
## luts are created if they do not exist
## reverse luts go from rpath -> aot
## QV 2021-02-03

def reverse_lut(sensor, lutdw=None, par = 'romix',
                       pct = (1,60), nbins = 20, override = False,
                       base_luts = ['ACOLITE-LUT-202101-MOD1', 'ACOLITE-LUT-202101-MOD2'],
                       rsky_base = 'ACOLITE-RSKY-202101-75W-{}ms', rsky_winds = [1,2,5,10]):
    import acolite as ac
    import numpy as np
    from netCDF4 import Dataset
    import scipy.interpolate
    import time, os

    if lutdw is None:
        rsrf = ac.config['data_dir']+'/RSR/{}.txt'.format(sensor)
        rsr, rsr_bands = ac.shared.rsr_read(rsrf)
        bands = [b for b in rsr_bands]
    else:
        bands = list(lutdw[lut]['rgi'].keys())

    lutdir = '{}-Reverse/{}'.format(ac.config['lut_dir'], sensor)
    if not os.path.exists(lutdir): os.makedirs(lutdir)

    revl = {}
    for lut in base_luts:
        rgi = {}
        for b in bands:
            slut = '{}-reverse-{}-{}-{}'.format(lut, sensor, par, b)
            lutnc = '{}/{}.nc'.format(lutdir, slut)

            if (not os.path.exists(lutnc)) or (override):
                if lutdw is None:
                    print('Importing source LUTs')
                    lutdw = ac.aerlut.import_luts(sensor=sensor, base_luts = base_luts,
                                                    add_rsky=True, add_rsky_winds=True,
                                                    rsky_base = rsky_base, rsky_winds = rsky_winds)
                    pid = lutdw[lut]['ipd'][par]
                    pressures, pids, raas, vzas, szas, winds, aots = lutdw[lut]['dim']

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
                                    ret = lutdw[lut]['rgi'][b]((pressure, pid, raa, vza, sza, wind, aots))
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
                var = nc.createVariable('lut',float,lut_dimensions)
                var[:] = luta
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
                    print(meta.keys())
                ## band specific interpolator
                rgi[b] = scipy.interpolate.RegularGridInterpolator([meta[k] for k in meta['lut_dimensions']],
                                                                 lutb,bounds_error=False, fill_value=None)
        revl[lut]={'rgi':rgi, 'minaot':minaot, 'maxaot':maxaot,
                   'model':int(lut[-1]), 'meta':meta}
    return(revl)
