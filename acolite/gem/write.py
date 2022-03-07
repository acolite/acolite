## write gem to L1R.nc file
## written by Quinten Vanhellemont, RBINS
## 2021-02-26
## modifications:

def write(gemfile, gem, verbosity=0):
    import os,json,bz2,dateutil.parser
    import numpy as np
    import acolite as ac

    ## global attributes
    keys = ['sensor', 'vza', 'sza', 'raa', 'isodate', 'pid']
    gatts = {k:gem[k] for k in keys}

    stime = dateutil.parser.parse(gatts['isodate'])
    doy = stime.strftime('%j')
    gatts['se_distance'] = ac.shared.distance_se(doy)

    oname = '{}_{}'.format(gatts['sensor'], stime.strftime('%Y_%m_%d_%H_%M_%S'))
    #if vname != '': oname+='_{}'.format(vname)

    #ofile = '{}/{}_gem_L1R.nc'.format(output, oname)
    ofile = '{}'.format(gemfile)
    if not os.path.exists(os.path.dirname(ofile)): os.makedirs(os.path.dirname(ofile))
    gatts['oname'] = oname
    gatts['ofile'] = ofile

    metakeys = list(gem['meta'].keys())
    metakeys.sort()
    for k in metakeys:
        if k in ['system:footprint']:
            continue
        else:
            gatts[k] = gem['meta'][k]

    ## add projection information
    pkeys = ['xrange', 'yrange', 'proj4_string', 'pixel_size', 'zone']
    for k in pkeys:
        if k in gem['dct']: gatts[k] = gem['dct'][k]

    ## write lon and lat
    ac.output.nc_write(ofile, 'lon', gem['data']['lon'], attributes=gatts, new=True)
    ac.output.nc_write(ofile, 'lat', gem['data']['lat'])

    for k in ['SAA', 'SZA', 'VAA', 'VZA', 'RAA']:
        if k in gem['data']:
            ac.output.nc_write(ofile, k.lower(), gem['data'][k])

    ## get band info
    rsrd = ac.shared.rsr_dict(sensor=gatts['sensor'])

    ## write bands
    for ib, b in enumerate(rsrd[gatts['sensor']]['rsr_bands']):
        btag = 'B{}'.format(b)
        ds = 'rhot_{}'.format(rsrd[gatts['sensor']]['wave_name'][b])
        ds_att  = {'wavelength':rsrd[gatts['sensor']]['wave_nm'][b]}
        ac.output.nc_write(ofile, ds, gem['data'][btag], attributes=gatts, dataset_attributes=ds_att)

    ## write thermal bands
    if gem['sensor'] in ['L5_TM', 'L7_ETM', 'L8_OLI', 'L9_OLI']:
        if  gem['sensor'] == 'L5_TM':
            thermal_bands = ['B6']
        if  gem['sensor'] == 'L7_ETM':
            thermal_bands = ['B6_VCID_1', 'B6_VCID_2']
        if  gem['sensor'] == 'L8_OLI':
            thermal_bands = ['B10', 'B11']
        if  gem['sensor'] == 'L9_OLI':
            thermal_bands = ['B10', 'B11']

        for btag in thermal_bands:
            if btag in gem['data']:
                ds = 'bt{}'.format(btag[1:].lower())
                ds_att  = {'K1': gem['meta']['K1_CONSTANT_BAND_{}'.format(btag[1:])],
                           'K2': gem['meta']['K2_CONSTANT_BAND_{}'.format(btag[1:])]}
                ac.output.nc_write(ofile, ds, gem['data'][btag])

    if verbosity > 1: print('Wrote {}'.format(gemfile))
    return(ofile)
