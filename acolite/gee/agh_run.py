## def agh_run
## runs ACOLITE/GEE hybrid/combined processing
## extracts data from GEE and does hybrid or local DSF
## written by Quinten Vanhellemont, RBINS
## 2022-04-13
## modifications: 2022-04-14 (QV) changed settings parsing, added option to add region name to output
##                2022-12-25 (QV) added SR option
##                2022-12-27 (QV) added offline TACT option
##                2023-06-21 (QV) use agh_old and added return_im keyword

def agh_run(settings={}, acolite_settings=None, rsrd = {}, lutd = {}, return_im=False, old_agh=True):
    import acolite as ac
    import os, time

    import ee
    #ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
    #ee.Authenticate() ## assume ee use is authenticated in current environment

    from acolite import gee

    ## get defaults
    s = ac.acolite.settings.read(ac.config['path']+'/config/gee_settings.txt')
    setg = ac.acolite.settings.parse(None, settings=s, merge=False)
    seti = ac.acolite.settings.parse(None, settings=settings, merge=False)
    for k in seti: setg[k] = seti[k]

    if (setg['add_region_name_output']) & (setg['region_name'] is not None):
        if setg['output'] is not None:
            setg['output'] = '{}/{}'.format(setg['output'], setg['region_name'])
        else:
            setg['output'] = '{}'.format(setg['region_name'])

    ## make our own limit if st_lat and st_lon are given with st_crop
    if (setg['st_lat'] is not None) & (setg['st_lon'] is not None) & (setg['st_crop']):
        setg['limit'] = ac.shared.region_box(None, setg['st_lon'], setg['st_lat'], \
                                             return_limit=True, box_size=setg['st_box'])
        print('Custom limit: ', setg['limit'])
        if setg['region_name_add_box'] & (setg['region_name'] is not None):
            setg['region_name'] += '_{:.0f}kmx{:.0f}km'.format(setg['st_box'],setg['st_box'])
            print('New region name: {}'.format(setg['region_name']))

    ## find images
    images, imColl = gee.find_scenes(setg['isodate_start'], isodate_end = setg['isodate_end'],
                                 day_range = setg['day_range'], sensors = setg['sensors'],
                                 surface_reflectance = setg['surface_reflectance'],
                                 filter_tiles = setg['filter_tiles'],
                                 limit = setg['limit'], st_lat = setg['st_lat'], st_lon = setg['st_lon'])

    print('Found images ', len(images), images)
    if return_im: return(images, imColl)

    ## run through images
    for image in images:
        skip = False
        if setg['filter_tiles'] is not None:
            if len(setg['filter_tiles'])>0:
                skip = True
                if any([tile in image[1] for tile in setg['filter_tiles']]): skip = False
        print(image, 'skip: ', skip)
        if skip: continue

        ## run AGH processing
        t0 = time.time()
        if old_agh:
            ret = gee.agh_old(image, imColl, rsrd=rsrd, lutd=lutd, settings=setg)
        else:
            ret = gee.agh(image, imColl, rsrd=rsrd, lutd=lutd, settings=setg)
        t1 = time.time()
        print('AGH processing finished in {:.1f} seconds'.format(t1-t0))

        ## run offline acolite
        if ('l1r_gee' in ret) & ('l2r_gee' not in ret) & (setg['run_offline_dsf'] | setg['run_offline_tact']):
            setu = ac.acolite.settings.parse(None, settings=acolite_settings, merge=False)
            setu['inputfile'] = ret['l1r_gee']
            if 'output' not in setu: setu['output'] = os.path.dirname(setu['inputfile'])
            print(setu)
            setu['atmospheric_correction'] = setg['run_offline_dsf']
            setu['tact_run'] = setg['run_offline_tact']
            ac.acolite.acolite_run(setu)
            t2 = time.time()
            print('Offline ACOLITE processing finished in {:.1f} seconds'.format(t2-t1))
