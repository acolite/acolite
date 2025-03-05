## def l1_convert
## converts Wyvern tiff file to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2025-03-04
## modifications: 2025-03-04 (QV) functionised from development code
##                2025-03-05 (QV) fixed issue with out of scene limit
##                                added wyvern_use_provided_f0 and wyvern_use_rsr_file

def l1_convert(inputfile, output = None, settings = None):
    import os, json
    import numpy as np
    import dateutil.parser
    import acolite as ac

    from osgeo import gdal
    gdal.UseExceptions()

    ## get run/user/sensor settings
    setu = ac.acolite.settings.merge(sensor = None, settings = settings)

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if setu['verbosity'] > 1: print('Starting conversion of {} scene{}'.format(nscenes, '' if nscenes==1 else 's'))

    ## list to store output files
    ofiles = []

    ## run through inputfiles
    for bi, bundle in enumerate(inputfile):

        ## get image and metadata file
        file, jf = ac.wyvern.bundle_test(bundle)
        meta, gatts = ac.wyvern.metadata_parse(jf)

        ## update settings
        setu = ac.acolite.settings.merge(sensor = gatts['sensor'], settings = settings)

        ## set output directory
        if output is None: output = setu['output']
        if output is None: output = os.path.dirname(file)

        ## get band info
        if setu['wyvern_use_rsr_file']: ## get RSR from file
            rsrd = ac.shared.rsr_dict(sensor = gatts['sensor'])
        else: ## get band info from gdal.Info
            info = gdal.Info(file).split('\n')
            band_data = {}
            new_band = False
            for il, line in enumerate(info):
                line = line.strip()
                if line.startswith('Band'):
                    new_band = True
                    cur_band = line.split(' ')[1]
                    band_data[cur_band] = {}
                    continue
                ## store band info
                if new_band:
                    sp = line.split('=')
                    if len(sp) == 2:
                        band_data[cur_band][sp[0]] = sp[1]
            ## add band info
            gatts['band_waves'] = [float(band_data[band]['wavelength']) for band in band_data]
            gatts['band_widths'] = [float(band_data[band]['FWHM']) for band in band_data]
            ## make rsr and bands dataset
            rsr = ac.shared.rsr_hyper(gatts['band_waves'], gatts['band_widths'], step=0.1)
            rsrd = ac.shared.rsr_dict(rsrd={gatts['sensor']:{'rsr':rsr}})

        ## read F0
        if setu['wyvern_use_provided_f0']:
            if setu['verbosity'] > 2: print('Using provided F0')
            f0d = {'{}'.format(bi+1): b['solar_illumination'] for bi, b in enumerate(meta['assets']['Cloud optimized GeoTiff']['eo:bands'])}
        else:
            if setu['verbosity'] > 2: print('Using F0 from solar_irradiance_reference={}'.format(setu['solar_irradiance_reference']))
            f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
            f0d = ac.shared.rsr_convolute_dict(f0['wave']/1000, f0['data'], rsrd[gatts['sensor']]['rsr'])

        ## make bands dataset
        bands = {}
        for bi, b in enumerate(rsrd[gatts['sensor']]['rsr']):
            cwave = rsrd[gatts['sensor']]['wave_nm'][b]
            swave = '{:.0f}'.format(cwave)
            bands[b]= {'wave':cwave, 'wavelength':cwave, 'wave_mu':cwave/1000.,
                       'wave_name':'{:.0f}'.format(cwave),
                       'rsr': rsrd[gatts['sensor']]['rsr'][b],'f0': f0d[b]}
            if 'band_widths' in gatts: bands[b]['width'] = gatts['band_widths'][bi]

        ## parse date
        time = dateutil.parser.parse(gatts['isodate'])
        doy = time.strftime('%j')
        se_distance = ac.shared.distance_se(doy)

        ## add to gatts
        gatts['doy'] = doy
        gatts['se_distance'] = se_distance
        gatts['acolite_file_type'] =  'L1R'

        ## output file name
        oname = '{}_{}'.format(gatts['sensor'], time.strftime('%Y_%m_%d_%H_%M_%S'))
        if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
        ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
        gatts['oname'] = oname
        gatts['ofile'] = ofile

        ## set up projection
        warp_to, dct_prj, sub = None, None, None
        try:
            ## get projection from image
            dct = ac.shared.projection_read(file)
        except:
            if setu['verbosity'] > 1: print('Could not determine image projection')
            dct = None

        ## find crop
        if (setu['limit'] is not None) and (dct is not None):
            dct_sub = ac.shared.projection_sub(dct, setu['limit'])
            if dct_sub['out_lon']:
                if setu['verbosity'] > 1: print('Longitude limits outside {}'.format(bundle))
                continue
            if dct_sub['out_lat']:
                if setu['verbosity'] > 1: print('Latitude limits outside {}'.format(bundle))
                continue
            sub = dct_sub['sub']

        if dct is not None:
            if sub is None:
                dct_prj = {k:dct[k] for k in dct}
            else:
                dct_prj = {k:dct_sub[k] for k in dct_sub}
                xyr = [min(dct_prj['xrange']),
                       min(dct_prj['yrange']),
                       max(dct_prj['xrange']),
                       max(dct_prj['yrange']),
                       dct_prj['proj4_string']]
                res_method = 'near'
                warp_to = (dct_prj['proj4_string'], xyr, dct_prj['pixel_size'][0],dct_prj['pixel_size'][1], res_method)

        ## for radiance conversion
        mus = np.cos(gatts['sza']*(np.pi/180))
        muv = np.cos(gatts['vza']*(np.pi/180))

        ## set up output gem
        gemo = ac.gem.gem(ofile, new = True)
        gemo.gatts = {k: gatts[k] for k in gatts}

        ## output geolocation
        if dct_prj is not None:
            if setu['verbosity'] > 1: print('Computing and writing lat/lon')
            ## offset half pixels to compute center pixel lat/lon
            dct_prj['xrange'] = dct_prj['xrange'][0]+dct_prj['pixel_size'][0]/2, dct_prj['xrange'][1]-dct_prj['pixel_size'][0]/2
            dct_prj['yrange'] = dct_prj['yrange'][0]+dct_prj['pixel_size'][1]/2, dct_prj['yrange'][1]-dct_prj['pixel_size'][1]/2
            ## compute lat/lon
            lon, lat = ac.shared.projection_geo(dct_prj, add_half_pixel = False)
            print(lat.shape)
            gemo.write('lat', lat)
            lat = None
            gemo.write('lon', lon)
            lon = None

        ## read data
        if setu['hyper_read_cube']:
            if setu['verbosity'] > 1: print('Reading Wyvern image cube')
            cube = ac.shared.read_band(file, sub = sub, warp_to = warp_to).astype(np.float32)
            print(cube.shape)

        ## run through bands
        for bi, b in enumerate(bands):
            if setu['verbosity'] > 2: print('Computing rhot_{} for {}'.format(bands[b]['wave_name'], gatts['oname']))
            ds_att = {k: bands[b][k] for k in bands[b] if k not in ['rsr']}

            wave = bands[b]['wavelength'] # gatts['band_waves'][bi]
            ds = 'rhot_{:.0f}'.format(wave)

            ## read data
            if setu['hyper_read_cube']:
                cdata_radiance = 1.0 * cube[bi, :, :]
            else:
                cdata_radiance = ac.shared.read_band(file, idx = bi+1, sub = sub, warp_to = warp_to).astype(np.float32)

            ## write toa radiance
            if setu['output_lt']:
                gemo.write('Lt_{}'.format(bands[b]['wave_name']), cdata_radiance, ds_att = ds_att)

            ## compute reflectance
            cdata = cdata_radiance * (np.pi * gatts['se_distance'] * gatts['se_distance']) / (bands[b]['f0'] * mus)
            cdata_radiance = None
            gemo.write('rhot_{}'.format(bands[b]['wave_name']), cdata, ds_att = ds_att)
            cdata = None
        cube = None
        gemo.close()

        ofiles.append(ofile)
    return(ofiles, setu)
