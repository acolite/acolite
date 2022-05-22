# def tact_gem
## runs tact on gem file
## written by Quinten Vanhellemont, RBINS
## 2021-02-28
## modifications: 2021-02-28 (QV) allow gem to be a dict
##                2022-03-22 (QV) added check whether thermal bands are in input gem

def tact_gem(gem, output_file = True,
             return_data = False,
             target_file = None,
             target_file_append = False,
             output = None,
             to_celcius = False,
             emissivity = 'water', # 'unity' / 'water' / 'eminet' / 'user'
             em_water = {'L5_TM':{'6': 0.9904561594813298},
                         'L7_ETM':{'6': 0.9909058218284508},
                         'L8_TIRS':{'10': 0.9926494385655781, '11': 0.9877384047023862},
                         'L9_TIRS':{'10': 0.9925715715215025, '11': 0.9876691339772199}},
             em_user = 1.0,
             sub = None,
             output_atmosphere = False,
             output_intermediate = False,
             copy_datasets = ['lon', 'lat'], verbosity=0):

    import os, datetime
    import numpy as np
    import acolite as ac

    ## read gem file if NetCDF
    if type(gem) is str:
        gemf = '{}'.format(gem)
        gem = ac.gem.read(gem, sub=sub)
    gemf = gem['gatts']['gemfile']

    ## check if we need to run tact
    run_tact = False
    for b in gem['gatts']['thermal_bands']:
        dsi = 'bt{}'.format(b)
        if dsi in gem['datasets']: run_tact = True
    if not run_tact:
        print('No thermal bands for {} (bands {}) in {}'.format(gem['gatts']['sensor'], ','.join(gem['gatts']['thermal_bands']), gemf))
        return()

    max_date = (datetime.datetime.now() - datetime.timedelta(days=90)).isoformat()
    if gem['gatts']['isodate'] > max_date:
        print('File too recent for TACT: after {}'.format(max_date))
        return()

    if 'nc_projection' in gem:
        nc_projection = gem['nc_projection']
    else:
        nc_projection = None

    if verbosity > 0: print('Running tact for {}'.format(gemf))

    ## detect sensor
    if ('thermal_sensor' not in gem['gatts']) or ('thermal_bands' not in gem['gatts']):
        if verbosity > 0: print('Processing of {} not supported'.format(gem['gatts']['sensor']))
        return()

    if target_file is None:
        #if 'output_name' in gem['gatts']:
        #    output_name = gem['gatts']['output_name']
        #elif 'oname' in gem['gatts']:
        #    output_name = gem['gatts']['oname']
        #else:
        #    output_name = os.path.basename(gemf).replace('.nc', '')
        output_name = os.path.basename(gemf).replace('.nc', '')
        output_name = output_name.replace('_L1R', '')
        output_name = output_name.replace(gem['gatts']['sensor'], gem['gatts']['thermal_sensor'])
        odir = output if output is not None else os.path.dirname(gemf)
        ofile = '{}/{}_ST.nc'.format(odir, output_name)
    else:
        ofile = '{}'.format(target_file)

    ## datasets to write
    output_datasets = []
    for ds in copy_datasets: output_datasets.append(ds)

    ## radiative transfer
    thd, simst, lonc, latc = ac.tact.tact_limit(gem['gatts']['isodate'],
                                                lon=gem['data']['lon'],
                                                lat=gem['data']['lat'],
                                                satsen=gem['gatts']['thermal_sensor'])
    for ds in thd:
        gem['data'][ds] = thd[ds]
        ## output atmosphere parameters
        if output_atmosphere: output_datasets += [ds]
    thd = None


    ## read bands and do thermal a/c
    for b in gem['gatts']['thermal_bands']:
        dsi = 'bt{}'.format(b)

        if dsi in gem['datasets']:
            btk = 'bt{}'.format(b)
            ltk = 'lt{}'.format(b)
            lsk = 'ls{}'.format(b)
            emk = 'em{}'.format(b)
            dso = 'st{}'.format(b)

            if output_intermediate: output_datasets += [btk, ltk, lsk, emk]
            output_datasets += [dso]

            #gd['data'][btk] = ac.shared.nc_data(ncf, dsi, sub=sub)
            #mask = gem['data'][btk].mask
            #gem['data'][btk] = gem['data'][btk].data
            #gem['data'][btk][mask] = np.nan

            bk = b.split('_')[0]
            if emissivity == 'unity':
                e = 1.0
            elif emissivity == 'water':
                e = em_water[gem['gatts']['thermal_sensor']][bk]
            elif emissivity == 'eminet':
                print('Emissivity from eminet to be implemented')
            elif emissivity == 'user':
                e = em_user
            else:
                continue

            ## shape emissivity to tile dimensions
            gem['data'][emk] = np.atleast_2d(e)
            if gem['data'][emk].shape == (1,1):
                gem['data'][emk] = np.tile(gem['data'][emk], gem['data']['lat'].shape)

            ## K constants
            k1 = float(gem['gatts']['K1_CONSTANT_BAND_{}'.format(b.upper())])
            k2 = float(gem['gatts']['K2_CONSTANT_BAND_{}'.format(b.upper())])

            ## compute lt from bt
            #bt = k2/(np.log(k1/lt)+1)
            #lt = k1/(np.exp(k2/bt)-1)
            gem['data'][ltk] = k1/(np.exp(k2/gem['data'][btk])-1)

            ## get surface radiance
            #ls = (((lt-thd['Lu{}'.format(bk)])/thd['tau{}'.format(bk)])-((1-e)*thd['Ld{}'.format(bk)]))/e
            gem['data'][lsk] = (((gem['data'][ltk]-gem['data']['Lu{}'.format(bk)])/gem['data']['tau{}'.format(bk)])-((1-gem['data'][emk])*gem['data']['Ld{}'.format(bk)]))/gem['data'][emk]

            ## convert to surface temperature
            #st = (k2/np.log((k1/ls)+1))
            gem['data'][dso] = (k2/np.log((k1/gem['data'][lsk])+1))
            if to_celcius: gem['data'][dso] += -273.15


    ## write output NetCDF
    if output_file:
        gem['gatts']['acolite_file_type'] = 'L2T'
        gem['gatts']['ofile'] = ofile
        gem['gatts']['sensor'] = gem['gatts']['thermal_sensor']

        new = True
        datasets_ofile = []
        if os.path.exists(ofile) & target_file_append:
            datasets_ofile = ac.shared.nc_datasets(ofile)
            new = False
        for ds in output_datasets:
            if ds in datasets_ofile: continue
            if ds not in gem['data']: continue
            ds_att = None
            if ds in gem['atts']: ds_att = gem['atts'][ds]
            ac.output.nc_write(ofile, ds, gem['data'][ds], new=new, nc_projection=nc_projection,
                               attributes=gem['gatts'], dataset_attributes=ds_att)
            if verbosity > 1: print('Wrote {} to {}'.format(ds, ofile))
            new=False
        if verbosity > 0: print('Wrote {}'.format(ofile))

    if return_data: return(gem)
    return(ofile)
