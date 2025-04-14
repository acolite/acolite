## def convert_gaac
## converts GAAC output bundle to netcdf file
## written by Quinten Vanhellemont, RBINS
## 2025-04-10
## modifications: 2025-04-14 (QV) add half pixel for lat lon computation

def convert_gaac(file, output = None, use_gaac_name = True):
    import acolite as ac
    import glob, os, re, dateutil.parser
    from osgeo import gdal
    gdal.UseExceptions()

    if os.path.isdir(file):
        dn = '{}'.format(file)
    else:
        dn = os.path.dirname(file)
    od = os.path.basename(dn)

    ## find scene files
    scene_files = glob.glob('{}/*'.format(dn))
    scene_files.sort()

    ## read tiff data and log
    log = []
    data = {}
    for sf in scene_files:
        if sf.endswith('.png'): continue ## skip png output
        if sf.endswith('rhor_rgb.tif'): continue ## skip rhor rgb output

        if sf.endswith('.txt'):
            with open(sf, 'r', encoding = 'utf-8') as f:
                for il, l in enumerate(f.readlines()):
                    log.append(l.strip())
        if sf.endswith('.tif'):
            ## get dataset names
            datasets = []
            with gdal.Open(sf) as ds:
                for bi in range(ds.RasterCount):
                    bds = ds.GetRasterBand(bi+1)
                    dataset_name = bds.GetDescription()
                    datasets.append(dataset_name)

            ## read file projection and determine lat/lon
            if sf.endswith('rhor.tif'):
                dct = ac.shared.projection_read(sf)
                lon, lat = ac.shared.projection_geo(dct, add_half_pixel = True)
                data['lon'] = lon
                data['lat'] = lat

            ## read tiff data
            meta, stack = ac.shared.read_band(sf, gdal_meta = True)
            if len(stack.shape) == 3:
                if len(datasets) == stack.shape[0]:
                    dataset_names = [ds for ds in datasets]
                    ## remove brackets from rhow datasets
                    for di, ds in enumerate(dataset_names):
                        if 'rhow' in ds:
                            dataset_names[di] = 'rhow_{}'.format(int(re.findall(r'\d+', ds)[0]))
                else:
                    dataset_names = ['band_{}'.format(i+1) for i in range(stack.shape[0])]
                for di, ds in enumerate(dataset_names):
                    data[ds] = stack[di, :, :]
            else:
                data[datasets[0]] = stack

    ## to do global attributes
    gatts = {}
    ## add rhor for spectrum viewer
    gatts['auto_grouping'] = 'rhot:rhorc:rhos:rhow:Rrs:Lt:Ed:rhor'

    ## extract info from log file
    if len(log) > 0:
        t0 = dateutil.parser.parse('T'.join(log[0].split()[0:2]).replace('/','-'))
        t1 = dateutil.parser.parse('T'.join(log[-1].split()[0:2]).replace('/','-'))
        gatts['processing_time'] = (t1-t0).total_seconds()

        par_keys = ['rhom', 'rhoa', 'trans_up', 'trans_down', 'trans_dir_up', 'trans_dir_down',
                    'Rrs', 'rhosky', 'rhog', 'rhoadj', 'q_factors', 'rhom_obs', 'Rrs_obs']
        oname = None

        for l in log:
            l = l.strip()
            if 'Output Dir:' in l:
                sp = l.split(':')
                gaac_name = sp[-1].strip()
                for v in gaac_name.split('_'):
                    if (len(v) >= 13) & ('T' in v):
                        dt = dateutil.parser.parse(v+'+00:00')
                        gatts['isodate'] = dt.isoformat()
                if use_gaac_name: oname = os.path.basename(gaac_name)

            if 'resolution:' in l:
                sp = l.split(':')
                gatts['resolution'] = float(sp[1])

            if 'sensor:' in l:
                sp = l.split(':')
                if sp[1] == 'S2_MSIA':
                    gatts['sensor'] = 'S2A_MSI'
                elif sp[1] == 'S2_MSIB':
                    gatts['sensor'] = 'S2B_MSI'
                elif sp[1] == 'S2_MSIC':
                    gatts['sensor'] = 'S2C_MSI'
                else:
                    print(sp)

            if 'Start GA optimization at pixel:' in l:
                sp = l.split(':')
                gatts['optimisation_location'] = sp[-1][1:sp[-1].find(']')]
            ## get additional keys
            for k in par_keys:
                if k in l:
                    sp = l.split()
                    sp = sp[-1].split(':')
                    gatts[k] = [float(v) for v in sp[1].split(',')]

    ## write nc
    odir = output if output is not None else os.path.dirname(dn)

    if oname is None:
        ofile = '{}/{}.nc'.format(odir, od)
    else:
        ofile = '{}/{}.nc'.format(odir, oname)
    gemo = ac.gem.gem(ofile, new = True)

    for ds in data:
        gemo.data_mem[ds] = data[ds]
        ## add wavelength for spectrum viewer
        if 'rho' in ds: gemo.data_att[ds] = {'wavelength': ds.split('_')[1]}
        gemo.write_ds(ds, clear = True)

    ## update attributes
    gemo.gatts = {k: gatts[k] for k in gatts}
    gemo.gatts_update()
    ## close file
    gemo.close()

    print('Wrote {}'.format(ofile))
    return(ofile)
