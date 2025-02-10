## def l1_convert
## converts IKONOS (old?) data to l1r NetCDF for acolite
## written by Quinten Vanhellemont, RBINS
## 2022-09-15
## modifications: 2022-11-21 (QV) added support for older/other metadata version from BELSPO
##                2023-07-12 (QV) removed netcdf_compression settings from nc_write call
##                2024-04-17 (QV) use new gem NetCDF handling
##                2025-02-02 (QV) removed percentiles
##                2025-02-04 (QV) improved settings handling
##                2025-02-10 (QV) cleaned up settings use, output naming

def l1_convert(inputfile, output = None, settings = None):

    import os, glob, dateutil.parser, datetime, time
    import numpy as np
    import scipy.ndimage
    import acolite as ac

    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}

    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    verbosity = setu['verbosity']

    ## parse inputfile
    if type(inputfile) != list:
        if type(inputfile) == str:
            inputfile = inputfile.split(',')
        else:
            inputfile = list(inputfile)
    nscenes = len(inputfile)
    if verbosity > 1: print('Starting conversion of {} scenes'.format(nscenes))

    new = True
    ofile = None
    ofiles = []

    for fi, bundle in enumerate(inputfile):
        sub = None
        t0 = time.time()

        if verbosity > 1: print('Importing metadata from {}'.format(bundle))
        files_dict = ac.ikonos.bundle_test(bundle)
        metadata = ac.ikonos.metadata_parse(files_dict['metadata'])

        ## try tags to use for retrieving components/sensor info
        try:
            ctag = 'Product Component Metadata'
            components = int(metadata[ctag]['Number of Components'][0])
            sensor = metadata['Source Image Metadata']['Sensor'][0].replace('-', '')
        except:
            ctag = 'Product Space Metadata'
            components = int(metadata[ctag]['Number of Image Components'][0])
            sensor = metadata['Product Order Metadata']['Sensor Name'][0].replace('-', '')

        ## get sensor specific defaults
        setd = ac.acolite.settings.parse(sensor)
        ## set sensor default if user has not specified the setting
        for k in setd:
            if k not in ac.settings['user']: setu[k] = setd[k]
        ## end set sensor specific defaults

        #if components != 1:
        #    print('Processing of IKONOS2 images with multiple components currently not supported.')
        #    continue
        verbosity = setu['verbosity']
        if output is None: output = setu['output']
        if output is None: output = os.path.dirname(bundle)

        ## read rsr
        rsrf = ac.path+'/data/RSR/{}.txt'.format(sensor)
        rsr, rsr_bands = ac.shared.rsr_read(rsrf)
        waves = np.arange(250, 2500)/1000
        waves_mu = ac.shared.rsr_convolute_dict(waves, waves, rsr)
        waves_names = {'{}'.format(b):'{:.0f}'.format(waves_mu[b]*1000) for b in waves_mu}

        gains = None
        if setu['gains']:
            if (len(setu['gains_toa']) == len(rsr_bands)) &\
               (len(setu['offsets_toa']) == len(rsr_bands)):
               gains = {}
               for bi, band in enumerate(rsr_bands):
                   gains[band] = {'gain': float(setu['gains_toa'][bi]),
                                'offset': float(setu['offsets_toa'][bi])}
            else:
                print('Use of gains requested, but provided number of gain ({}) or offset ({}) values does not match number of bands in RSR ({})'.format(len(setu['gains_toa']), len(setu['offsets_toa']), len(rsr_bands)))
                print('Provide gains in band order: {}'.format(','.join(rsr_bands)))

        ## get F0
        f0 = ac.shared.f0_get(f0_dataset=setu['solar_irradiance_reference'])
        f0_b = ac.shared.rsr_convolute_dict(np.asarray(f0['wave'])/1000, np.asarray(f0['data'])*10, rsr)

        ## global scene dimensions from metadata
        #global_dims = [int(metadata['Product Space Metadata']['Rows'].split(' ')[0]),
        #               int(metadata['Product Space Metadata']['Columns'].split(' ')[0])]

        ## write results to output file
        for ci in range(components):
            for t in ['Component ID', 'Tile ID']:
                print(ctag, t)
                try:
                    comp = metadata[ctag][t][ci]
                    break
                except:
                    pass

            dtime = dateutil.parser.parse(metadata['Source Image Metadata']['Acquisition Date/Time'][ci])
            isodate = dtime.isoformat()
            doy = dtime.strftime('%j')
            se_distance = ac.shared.distance_se(doy)

            ## IKONOS DN radiance calibration info
            # IKONOS Radiometric Calibration Coefficients [mW/(cm2*sr*DN)]
            # move to external file?
            if metadata['Product Order Metadata']['Bits per Pixel per Band'][0] == '11 bits per pixel':
                if isodate <= '2001-02-22':
                    cal = {'pan': 161, 'blu': 633, 'grn': 649, 'red': 840, 'nir': 746}
                elif isodate > '2001-02-22':
                    cal = {'pan': 161, 'blu': 728, 'grn': 727, 'red': 949, 'nir': 843}
            elif metadata['Product Order Metadata']['Bits per Pixel per Band'][0] == '8 bits per pixel':
                if isodate <= '2001-02-22':
                    cal = {'pan': 161, 'blu': 79, 'grn': 81, 'red': 105, 'nir': 93}
                elif isodate > '2001-02-22':
                    cal = {'pan': 161, 'blu': 91, 'grn': 91, 'red': 119, 'nir': 105}

            # IKONOS band widths
            bandwidths = {'pan': 403, 'blu': 71.3, 'grn': 88.6, 'red': 65.8, 'nir': 95.4}

            ## get observation geometry
            vaa = float(metadata['Source Image Metadata']['Nominal Collection Azimuth'][ci].split(' ')[0])
            vea = float(metadata['Source Image Metadata']['Nominal Collection Elevation'][ci].split(' ')[0])
            vza = 90 - vea
            saa = float(metadata['Source Image Metadata']['Sun Angle Azimuth'][ci].split(' ')[0])
            sea = float(metadata['Source Image Metadata']['Sun Angle Elevation'][ci].split(' ')[0])
            sza = 90 - sea
            raa = abs(saa - vaa)
            while raa >= 180.: raa = np.abs(raa-360)

            gatts = {'sensor': sensor, 'satellite': sensor,'isodate': isodate,
                         'sza': sza, 'vza': vza, 'saa': saa, 'vaa': vaa,
                         'raa': raa, 'se_distance': se_distance,
                         'mus': np.cos(sza*(np.pi/180.)), 'acolite_file_type': 'L1R'}
            ## add band info to gatts
            for b in rsr_bands:
                gatts['{}_wave'.format(b)] = waves_mu[b]*1000
                gatts['{}_name'.format(b)] = waves_names[b]
                gatts['{}_f0'.format(b)] = f0_b[b]

            stime = dateutil.parser.parse(gatts['isodate'])

            oname = '{}_{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'), comp)
            if setu['region_name'] != '': oname+='_{}'.format(setu['region_name'])
            ofile = '{}/{}_{}.nc'.format(output, oname, gatts['acolite_file_type'])
            pofile = '{}/{}_{}_pan.nc'.format(output, oname, gatts['acolite_file_type'])
            gatts['oname'] = oname
            gatts['ofile'] = ofile
            gatts['ofile_pan'] = pofile

            new=True
            new_pan = True

            band_names = ['Blue', 'Green', 'Red', 'NIR', 'PAN']
            band_keys = ['blu', 'grn', 'red', 'nir', 'pan']

            for bi, band in enumerate(band_names):
                bk = band_keys[bi]
                if bk not in files_dict: continue
                band_file = None
                for bf in files_dict[bk]:
                    if comp in os.path.basename(bf):
                        band_file = '{}'.format(bf)
                        break

                if band_file == None:
                    print('Band file for band {} component {} not found'.format(band, comp))
                    continue

                ## read band image projection
                ## for some files only the blu has the projection info?
                try:
                    dct_cur = ac.shared.projection_read(band_file)
                except:
                    print('Could not determine projection for band {} {}'.format(band, band_file))

                    pass

                if bi == 0:
                    dct = {k: dct_cur[k] for k in dct_cur}
                    nc_projection = ac.shared.projection_netcdf(dct)
                    global_dims = len(nc_projection['y']['data']), len(nc_projection['x']['data'])

                if (bk != 'pan') & (dct != dct_cur):
                    print('Warning, band {} has a different projection dct'.format(band))

                if new:
                    ## output file - one per component
                    gemo = ac.gem.gem(ofile, new = True)
                    gemo.gatts = {k: gatts[k] for k in gatts}
                    gemo.nc_projection = nc_projection
                    new = False

                ## write lat/lon
                if (setu['output_geolocation']) & (new):
                    if verbosity > 1: print('{} - Writing lat/lon'.format(datetime.datetime.now().isoformat()[0:19]))
                    if dct is not None: ## compute from projection info
                        lon, lat = ac.shared.projection_geo(dct, add_half_pixel=False)
                        gemo.write('lon', lon)
                        lon = None
                        gemo.write('lat', lat)
                        lat = None

                ## read band data
                if verbosity > 1: print('{} - Reading band {} from {}'.format(datetime.datetime.now().isoformat()[0:19],
                                                                              band, band_file))
                data = ac.shared.read_band(band_file)
                nodata = data == np.uint16(0)

                ## radiance
                cf = (10**4) / (cal[bk] * bandwidths[bk])
                data = data.astype(np.float32) * cf
                data[nodata] = np.nan

                ## apply gains
                if (gains != None) & (setu['gains_parameter'] == 'radiance'):
                    print('Applying gain {} and offset {} to TOA radiance for band {}'.format(gains[band]['gain'], gains[band]['offset'], band))
                    data = gains[band]['gain'] * data + gains[band]['offset']

                ## reflectance
                data *= (np.pi * gatts['se_distance']**2) / (f0_b[band]/10. * gatts['mus'])

                ## set up dataset attributes
                ds = 'rhot_{}'.format(waves_names[band])
                ds_att = {'wavelength': waves_mu[band]*1000, 'band_name': band, 'f0': f0_b[band]/10.}
                if gains != None:
                    ds_att['gain'] = gains[band]['gain']
                    ds_att['offset'] = gains[band]['offset']
                    ds_att['gains_parameter'] = setu['gains_parameter']

                ## save PAN separately and downsample
                if band == 'PAN':
                    dct_pan = {k: dct_cur[k] for k in dct_cur}
                    nc_projection_pan = ac.shared.projection_netcdf(dct_pan, add_half_pixel=False)
                    if new_pan:
                        gemop = ac.gem.gem(pofile, new = True)
                        gemop.gatts = {k: gatts[k] for k in gatts}
                        gemop.nc_projection = nc_projection_pan
                        new_pan = False
                        global_dims_pan = len(nc_projection_pan['y']['data']), len(nc_projection_pan['x']['data'])
                    ## write to netcdf file
                    gemop.write(ds, data, ds_att = ds_att, replace_nan = True)

                    ## mask data before zooming
                    dmin = np.nanmin(data)
                    data[np.isnan(data)] = 0
                    data = scipy.ndimage.zoom(data, 0.25, order=1)
                    data[data<dmin] = np.nan
                    data[data==dmin] = np.nan

                ## write to netcdf file
                if verbosity > 1: print('{} - Converting bands: Writing {} ({})'.format(datetime.datetime.now().isoformat()[0:19], ds, data.shape))
                gemo.write(ds, data, ds_att = ds_att)
                if verbosity > 1: print('{} - Converting bands: Wrote {} ({})'.format(datetime.datetime.now().isoformat()[0:19], ds, data.shape))
                data = None

            if not os.path.exists(ofile): continue
            if verbosity > 1:
                print('Conversion took {:.1f} seconds'.format(time.time()-t0))
                print('Created {}'.format(ofile))
            gemo.close()
            try:
                gemop.close()
            except:
                pass

            if ofile not in ofiles: ofiles.append(ofile)

    return(ofiles, setu)
