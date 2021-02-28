# def tact_gem
## runs tact on gem file
## written by Quinten Vanhellemont, RBINS
## 2021-02-28
## modifications:

def tact_gem(ncf, output_file = True,
             return_data = False,
             target_file = None,
             target_file_append = False,
             output = None,
             to_celcius = False,
             emissivity = 'water', # 'unity' / 'water' / 'eminet' / 'user'
             em_water = {'L5_TM':{'6': 0.9904561594813298},
                         'L7_ETM':{'6': 0.9909058218284508},
                         'L8_TIRS':{'10': 0.9926494385655781, '11': 0.9877384047023862}},
             em_user = 1.0,
             sub = None,
             output_atmosphere = False,
             output_intermediate = False,
             copy_datasets = ['lon', 'lat'], verbosity=0):

    import os
    import numpy as np
    import scipy.ndimage
    import acolite as ac

    if verbosity > 0: print('Running tact for {}'.format(ncf))

    ## get datasets and attributes from NetCDF
    datasets = ac.shared.nc_datasets(ncf)
    gatts = ac.shared.nc_gatts(ncf)

    ## detect sensor
    if gatts['sensor'] == 'L8_OLI':
        tact_sen = 'L8_TIRS'
        th_bands = ['10', '11']
    elif gatts['sensor'] == 'L5_TM':
        tact_sen = gatts['sensor']
        th_bands = ['6']
    elif gatts['sensor'] == 'L7_ETM':
        tact_sen = gatts['sensor']
        th_bands = ['6_vcid_1', '6_vcid_2']
    else:
        if verbosity > 0: print('Processing of {} not supported'.format(gatts['sensor']))
        return()

    if target_file is None:
        if 'output_name' in gatts:
            output_name = gatts['output_name']
        else:
            output_name = os.path.basename(ncf).replace('.nc', '')
        output_name = output_name.replace('_L1R', '')
        odir = output if output is not None else os.path.dirname(ncf)
        ofile = '{}/{}_ST.nc'.format(odir, output_name)
    else:
        ofile = '{}'.format(target_file)

    ## datasets to write
    output_datasets = []
    for ds in copy_datasets: output_datasets.append(ds)

    ## empty gem dict
    gd = {'data':{}, 'atts':{}}

    ## lat lon are needed to get profiles for tact
    for ds in ['lat', 'lon']:
        if ds in datasets:
            gd['data'][ds] = ac.shared.nc_data(ncf, ds, sub=sub)
            gd['atts'][ds] = None

    ## radiative transfer
    thd, simst, lonc, latc = ac.tact.tact_limit(gatts['isodate'],
                                                lon=gd['data']['lon'],
                                                lat=gd['data']['lat'],
                                                satsen=tact_sen)
    for ds in thd:
        gd['data'][ds] = thd[ds]
        ## output atmosphere parameters
        if output_atmosphere: output_datasets += [ds]
    thd = None


    ## read bands and do thermal a/c
    for b in th_bands:
        ## old (gee) gems have BT_B10
        dsi = 'BT_B{}'.format(b.upper())

        ## new gems have bt10
        if dsi not in datasets: dsi = 'bt{}'.format(b)

        if dsi in datasets:
            btk = 'bt{}'.format(b)
            ltk = 'lt{}'.format(b)
            lsk = 'ls{}'.format(b)
            emk = 'em{}'.format(b)
            dso = 'st{}'.format(b)

            if output_intermediate: output_datasets += [btk, ltk, lsk, emk]
            output_datasets += [dso]

            gd['data'][btk] = ac.shared.nc_data(ncf, dsi, sub=sub)
            mask = gd['data'][btk].mask
            gd['data'][btk] = gd['data'][btk].data
            gd['data'][btk][mask] = np.nan

            bk = b.split('_')[0]
            if emissivity == 'unity':
                e = 1.0
            elif emissivity == 'water':
                e = em_water[tact_sen][bk]
            elif emissivity == 'eminet':
                print('Emissivity from eminet to be implemented')
            elif emissivity == 'user':
                e = em_user
            else:
                continue

            ## shape emissivity to tile dimensions
            gd['data'][emk] = np.atleast_2d(e)
            if gd['data'][emk].shape == (1,1):
                gd['data'][emk] = np.tile(gd['data'][emk], gd['data']['lat'].shape)

            ## K constants
            k1 = float(gatts['K1_CONSTANT_BAND_{}'.format(b.upper())])
            k2 = float(gatts['K2_CONSTANT_BAND_{}'.format(b.upper())])

            ## compute lt from bt
            #bt = k2/(np.log(k1/lt)+1)
            #lt = k1/(np.exp(k2/bt)-1)
            gd['data'][ltk] = k1/(np.exp(k2/gd['data'][btk])-1)

            ## get surface radiance
            #ls = (((lt-thd['Lu{}'.format(bk)])/thd['tau{}'.format(bk)])-((1-e)*thd['Ld{}'.format(bk)]))/e
            gd['data'][lsk] = (((gd['data'][ltk]-gd['data']['Lu{}'.format(bk)])/gd['data']['tau{}'.format(bk)])-((1-gd['data'][emk])*gd['data']['Ld{}'.format(bk)]))/gd['data'][emk]

            ## convert to surface temperature
            #st = (k2/np.log((k1/ls)+1))
            gd['data'][dso] = (k2/np.log((k1/gd['data'][lsk])+1))
            if to_celcius: gd['data'][dso] += -273.15


    ## write output NetCDF
    if output_file:
        new = True
        datasets_ofile = []
        if os.path.exists(ofile) & target_file_append:
            datasets_ofile = ac.shared.nc_datasets(ofile)
            new = False
        for ds in output_datasets:
            if ds in datasets_ofile: continue
            if ds not in gd['data']: continue
            ds_att = None
            if ds in gd['atts']: ds_att = gd['atts'][ds]
            ac.output.nc_write(ofile, ds, gd['data'][ds], new=new, attributes=gatts, dataset_attributes=ds_att)
            if verbosity > 1: print('Wrote {} to {}'.format(ds, ofile))
            new=False
        if verbosity > 0: print('Wrote {}'.format(ofile))

    if return_data: return(gd)
    return(ofile)
