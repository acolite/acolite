# def tact_gem
## runs tact on gem file
## written by Quinten Vanhellemont, RBINS
## 2021-02-28
## modifications: 2021-02-28 (QV) allow gem to be a dict
##                2022-03-22 (QV) added check whether thermal bands are in input gem
##                2022-07-31 (QV) skip loading of datasets that are not required
##                2022-08-02 (QV) added source keyword
##                2022-08-03 (QV) added external emissivity files
##                2024-03-14 (QV) update settings handling
##                2024-04-17 (QV) use new gem NetCDF handling
##                2025-02-04 (QV) updated settings parsing
##                2025-02-10 (QV) cleaned up settings use

def tact_gem(gem, output_file = True,
             output = None,
             return_data = False,
             target_file = None,
             target_file_append = False,
             to_celcius = False,
             sub = None,
             copy_datasets = ['lon', 'lat'],
             settings = None, verbosity=0):

    import os, datetime, json
    import numpy as np
    import acolite as ac

    ## read gem file if NetCDF
    if type(gem) is str:
        gem = ac.gem.gem(gem)
    gemf = gem.file

    ## find datasets we need
    skip_datasets = []
    for ds in gem.datasets:
        if ds in ['lat', 'lon']: continue
        if ds[0:2] == 'bt': continue
        if ds[0:2] == 'lt': continue
        skip_datasets.append(ds)

    ## read in datasets from NetCDF
    dct = {'data': {}, 'atts': {}}
    for ds in gem.datasets:
        if ds in skip_datasets: continue
        if 'projection_key' in gem.gatts:
                if ds in ['x', 'y', gem.gatts['projection_key']]: continue
        d_, a_ = gem.data(ds, sub = sub, attributes = True)
        dct['data'][ds] = d_
        dct['atts'][ds] = a_
        print('Read {} {}'.format(ds, d_.shape))
        d_, a_= None, None

    ## combine default and user defined settings
    ## get run settings
    setu = {k: ac.settings['run'][k] for k in ac.settings['run']}
    ## get sensor specific defaults
    setd = ac.acolite.settings.parse(gem.gatts['sensor'])
    ## set sensor default if user has not specified the setting
    for k in setd:
        if k not in ac.settings['user']: setu[k] = setd[k]
    ## end set sensor specific defaults
    ## additional run settings
    if settings is not None:
        settings = ac.acolite.settings.parse(settings)
        for k in settings: setu[k] = settings[k]
    ## end additional run settings

    ## get verbosity from run settings
    verbosity = setu['verbosity']
    output = setu['output']

    ## detect sensor
    if ('thermal_sensor' not in gem.gatts) or ('thermal_bands' not in gem.gatts):
        if verbosity > 0: print('TACT Processing of {} not supported'.format(gem.gatts['sensor']))
        return

    ## check if we need to run tact
    run_tact = False
    for b in gem.gatts['thermal_bands']:
        dsi = 'bt{}'.format(b)
        if dsi in gem.datasets: run_tact = True
    if not run_tact:
        print('No thermal bands for {} (bands {}) in {}'.format(gem.gatts['sensor'], ','.join(gem.gatts['thermal_bands']), gemf))
        return

    ## check blackfill
    if setu['blackfill_skip']:
        for b in gem.gatts['thermal_bands']:
            if 'bt{}'.format(b) in dct['data']:
                band_data = 1.0*dct['data']['bt{}'.format(b)]
                break
        npx = band_data.shape[0] * band_data.shape[1]
        #nbf = npx - len(np.where(np.isfinite(band_data))[0])
        nbf = npx - len(np.where(np.isfinite(band_data)*(band_data>0))[0])
        band_data = None
        if (nbf/npx) >= float(setu['blackfill_max']):
            if verbosity>0: print('Skipping scene as crop is {:.0f}% blackfill'.format(100*nbf/npx))
            return

    if setu['tact_profile_source'] == 'era5':
        max_date = (datetime.datetime.now() - datetime.timedelta(days=92)).isoformat()
        if gem.gatts['isodate'] > max_date:
            print('File too recent for TACT with {} profiles: after {}'.format(setu['tact_profile_source'], max_date))
            print('Run with tact_profile_source=gdas1 or tact_profile_source=ncep.reanalysis2 for NRT processing')
            return

    ## load emissivity data
    emissivity_file = setu['tact_emissivity_file']
    if emissivity_file is not None:
        if not os.path.exists(emissivity_file):
            print('Could not file {}'.format(emissivity_file))
            emissivity_file = None
    if (emissivity_file is None) & (setu['tact_emissivity'] not in ['ged', 'eminet', 'ndvi']):
        emissivity_file = '{}/{}/emissivity_{}.json'.format(ac.config['data_dir'], 'TACT', setu['tact_emissivity'])
        if not os.path.exists(emissivity_file):
            print('Could not file {}'.format(emissivity_file))
            emissivity_file = None
    if emissivity_file is not None:
        em = json.load(open(emissivity_file, 'r'))
        print('Loaded emissivity file {}'.format(emissivity_file))
    else:
        em = None

    if verbosity > 0: print('Running tact for {}'.format(gemf))

    if target_file is None:
        output_name = os.path.basename(gemf).replace('.nc', '')
        output_name = output_name.replace('_L1R', '')
        output_name = output_name.replace(gem.gatts['sensor'], gem.gatts['thermal_sensor'])
        odir = output if output is not None else os.path.dirname(gemf)
        ofile = '{}/{}_ST.nc'.format(odir, output_name)
    else:
        ofile = '{}'.format(target_file)

    print('Running ACOLITE/TACT for {}'.format(gemf))

    ## read lon/lat
    lon = gem.data('lon')
    lat = gem.data('lat')

    ## datasets to write
    output_datasets = []
    for ds in copy_datasets: output_datasets.append(ds)

    ## radiative transfer
    thd, simst, lonc, latc = ac.tact.tact_limit(gem.gatts['isodate'],
                                                lon = lon, lat = lat,
                                                satsen=gem.gatts['thermal_sensor'],
                                                wave_range = setu['tact_range'],
                                                reptran = setu['tact_reptran'],
                                                source = setu['tact_profile_source'])
    for ds in thd:
        dct['data'][ds] = thd[ds]
        ## output atmosphere parameters
        if setu['tact_output_atmosphere']: output_datasets += [ds]
    thd = None


    ## read bands and do thermal a/c
    em_ged, em_eminet, em_ndvi = None, None, None
    for b in gem.gatts['thermal_bands']:
        dsi = 'bt{}'.format(b)

        if dsi in gem.datasets:
            btk = 'bt{}'.format(b)
            ltk = 'lt{}'.format(b)
            lsk = 'ls{}'.format(b)
            emk = 'em{}'.format(b)
            dso = 'st{}'.format(b)

            if setu['tact_output_intermediate']: output_datasets += [btk, ltk, lsk, emk]
            output_datasets += [dso]

            bk = b.split('_')[0]
            e = None
            if setu['tact_emissivity'] == 'ged':
                if em_ged is None:
                    ## determine bands
                    if gem.gatts['thermal_sensor'] in ['L8_TIRS', 'L9_TIRS']:
                        bands = [13, 14]
                        bkeys = {'10':0, '11':1}
                    elif gem.gatts['thermal_sensor'] in ['L5_TM', 'L7_ETM']:
                        bands = [13]
                        bkeys = {'6':0}
                    elif gem.gatts['thermal_sensor'] in ['ISS_ECOSTRESS']:
                        bands = [10, 11, 12, 13, 14]
                        bkeys = {'1':0, '2':1, '3':2, '4':3, '5':4}
                    ## load GED emissivity
                    em_ged = ac.ged.ged_lonlat(lon, lat, bands=bands, fill = setu['ged_fill'])
                if em_ged is None:
                    print('Could not extract GED emissivity.')
                else:
                    if len(em_ged.shape) == 3:
                        e = em_ged[:,:,bkeys[b]]
                    else:
                        e = em_ged * 1.0
            if setu['tact_emissivity'] == 'eminet':
                if em_eminet is None:
                    em_eminet = ac.tact.tact_eminet(gemf, water_fill = setu['eminet_water_fill'],
                                                       water_threshold = setu['eminet_water_threshold'],
                                                       model_version = setu['eminet_model_version'],
                                                       netname = setu['eminet_netname'],
                                                       fill = setu['eminet_fill'],
                                                       fill_dilate = setu['eminet_fill_dilate'])
                if em_eminet is None:
                    print('Could not get EMINET emissivity.')
                else:
                    e = em_eminet[bk]

            if setu['tact_emissivity'] == 'ndvi':
                if em_ndvi is None:
                    em_ndvi = ac.tact.ndvi_emissivity(gemf, ndvi_toa = setu['tact_emissivity_ndvi_toa'])
                if em_ndvi is None:
                    print('Could not get NDVI derived emissivity.')
                else:
                    e = em_ndvi[bk]

            if (e is None) & (em is not None):
                e = em[gem.gatts['thermal_sensor']][bk]
                #print(e, gem.gatts['thermal_sensor'], bk)

            if e is None:
                print('Emissivity for {} {} not configured.'.format(gem.gatts['thermal_sensor'], bk))
                print('Ls and ST will not be computed for {} {}.'.format(gem.gatts['thermal_sensor'], bk))
                continue

            ## shape emissivity to tile dimensions
            dct['data'][emk] = np.atleast_2d(e)
            if dct['data'][emk].shape == (1,1):
                dct['data'][emk] = np.tile(dct['data'][emk], dct['data']['lat'].shape)

            ## K constants
            k1n = 'K1_CONSTANT_BAND_{}'.format(b.upper())
            k2n = 'K2_CONSTANT_BAND_{}'.format(b.upper())
            if k1n in gem.gatts:
                k1 = float(gem.gatts[k1n])
            else:
                k1 = dct['atts'][btk][k1n]
            if k2n in gem.gatts:
                k2 = float(gem.gatts[k2n])
            else:
                k2 = dct['atts'][btk][k2n]

            if ltk not in dct['data']:
                ## compute lt from bt
                #bt = k2/(np.log(k1/lt)+1)
                #lt = k1/(np.exp(k2/bt)-1)
                dct['data'][ltk] = k1/(np.exp(k2/dct['data'][btk])-1)

            ## get surface radiance
            #ls = (((lt-thd['Lu{}'.format(bk)])/thd['tau{}'.format(bk)])-((1-e)*thd['Ld{}'.format(bk)]))/e
            dct['data'][lsk] = (((dct['data'][ltk]-dct['data']['Lu{}'.format(bk)])/dct['data']['tau{}'.format(bk)])-((1-dct['data'][emk])*dct['data']['Ld{}'.format(bk)]))/dct['data'][emk]

            ## convert to surface temperature
            #st = (k2/np.log((k1/ls)+1))
            dct['data'][dso] = (k2/np.log((k1/dct['data'][lsk])+1))
            if to_celcius: dct['data'][dso] += -273.15

            dct['atts'][dso] = {}
            dct['atts'][dso]['units'] = 'K'
            if btk in dct['atts']:
                if 'wavelength' in dct['atts'][btk]:
                    dct['atts'][dso]['wavelength'] = dct['atts'][btk]['wavelength']
            #if btk in dct['atts']:
            #    dct['atts'][dso] = {k:dct['atts'][btk][k] for k in dct['atts'][btk]}
            dct['atts'][dso][k1n] = k1
            dct['atts'][dso][k2n] = k2


    ## write output NetCDF
    if output_file:
        gatts = {k: gem.gatts[k] for k in gem.gatts}
        nc_projection = gem.nc_projection

        ## update gatts
        gatts['acolite_file_type'] = 'L2T'
        gatts['ofile'] = ofile
        gatts['sensor'] = gatts['thermal_sensor']

        new = True
        datasets_ofile = []
        if os.path.exists(ofile) & target_file_append:
            datasets_ofile = ac.shared.nc_datasets(ofile)
            new = False

        ## create new file
        if new:
            gemo = ac.gem.gem(ofile, new = True)
            gemo.gatts = {k: gatts[k] for k in gatts}
            gemo.nc_projection = nc_projection
            new = False
        else:
            gemo = ac.gem.gem(ofile)

        for ds in output_datasets:
            if ds in datasets_ofile: continue
            if ds not in dct['data']: continue
            ds_att = None
            if ds in dct['atts']: ds_att = dct['atts'][ds]
            gemo.write(ds, dct['data'][ds], ds_att = ds_att)
            if verbosity > 1: print('Wrote {} to {}'.format(ds, ofile))
        if verbosity > 0: print('Wrote {}'.format(ofile))
        gemo.close()
        gemo = None

    ## close gem
    gem.close()
    gem = None

    if return_data: return(dct)
    dct = None
    return(ofile)
