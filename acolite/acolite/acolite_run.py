## def acolite_run
## runs acolite processing for given settings file/dict, new version
## written by Quinten Vanhellemont, RBINS
## 2021-04-01
## modifications: 2021-04-14 (QV) added output to settings if not configured
##                2021-04-15 (QV) test/parse input files
##                2022-03-04 (QV) moved inputfile testing to inputfile_test
##                2022-07-25 (QV) avoid deleting original inputfiles
##                2022-09-19 (QV) printout platform info
##                2023-02-06 (QV) added WKT polygon output
##                2024-03-14 (QV) update settings handling
##                2024-03-28 (QV) added station limit creation
##                2024-13-17 (QV) update for atmospheric_correction_method
##                2025-01-20 (QV) updated RGB output logic
##                2025-01-30 (QV) moved polygon limit and limit buffer extension
##                2025-02-02 (QV) moved version printout up, fixed limit_buffer when limit is None

def acolite_run(settings, inputfile=None, output=None):
    import glob, datetime, os, shutil, copy
    import acolite as ac

    ## time of processing start
    time_start = datetime.datetime.now()

    ## reset run settings to defaults
    ac.settings['run'] = {k:ac.settings['defaults'][k] for k in ac.settings['defaults']}
    if 'runid' not in ac.settings['run']: ac.settings['run']['runid'] = time_start.strftime('%Y%m%d_%H%M%S')

    ## get user settings
    ## these are updated with sensor specific settings in acolite_l2r/l2w/tact
    ac.settings['user'] = ac.acolite.settings.parse(None, settings = settings, merge=False)
    if 'output' not in ac.settings['user']: ac.settings['run']['output'] = os.getcwd()

    ## set settings from launch_acolite
    if inputfile is not None: ac.settings['user']['inputfile'] = inputfile
    if output is not None: ac.settings['user']['output'] = output

    ## workaround for outputting rhorc and bt
    if 'l2w_parameters' in ac.settings['user']:
        if ac.settings['user']['l2w_parameters'] is not None:
            for par in ac.settings['user']['l2w_parameters']:
                if par.startswith('rhorc'): ac.settings['user']['output_rhorc'] = True
                if par.startswith('bt'): ac.settings['user']['output_bt'] = True
                if par.startswith('Ed'): ac.settings['user']['output_ed'] = True

    ## update run settings
    for k in ac.settings['user']: ac.settings['run'][k] = ac.settings['user'][k]
    if 'verbosity' in ac.settings['run']: ac.config['verbosity'] = int(ac.settings['run']['verbosity'])

    ## new path is only set if ACOLITE needs to make new directories
    ## and is only used if ACOLITE is asked to delete the output directory (don't use this feature!)
    new_path = None
    if 'new_path' in ac.settings['run']: new_path = '{}'.format(ac.settings['run']['new_path'])

    ## log file for l1r generation
    log_file = '{}/acolite_run_{}_log_file.txt'.format(ac.settings['run']['output'], ac.settings['run']['runid'])
    log = ac.acolite.logging.LogTee(log_file)

    print('Running ACOLITE processing - {}'.format(ac.version))
    print('Python - {} - {}'.format(ac.python['platform'], ac.python['version']).replace('\n', ''))
    print('Platform - {} {} - {} - {}'.format(ac.system['sysname'], ac.system['release'], ac.system['machine'], ac.system['version']).replace('\n', ''))
    print('Run ID - {}'.format(ac.settings['run']['runid']))

    ## check if we have anything to do
    if 'inputfile' not in ac.settings['run']:
        print('Nothing to do. Did you provide a settings file or inputfile? Exiting.')
        return()
    else:
        nscenes = len(ac.settings['run']['inputfile'])

    if ac.settings['run']['polylakes']:
        ac.settings['run']['polygon'] = ac.shared.polylakes(ac.settings['run']['polylakes_database'])
        ac.settings['run']['polygon_limit'] = False
        ac.settings['run']['polygon_clip'] = True

    ## check if polygon is a file or WKT
    ## if WKT make new polygon file in output directory
    ## new and old polygon inputs tracked here in user settings since WKT is provided by the user
    if 'polygon' in ac.settings['run']:
        if ac.settings['run']['polygon'] is not None:
            ## is the given polygon a wkt?
            if not os.path.exists(ac.settings['run']['polygon']):
                try:
                    polygon_new = ac.shared.polygon_from_wkt(ac.settings['run']['polygon'],
                                     file = '{}/polygon_{}.json'.format(ac.settings['run']['output'], ac.settings['run']['runid']))
                    ac.settings['run']['polygon_old'] = '{}'.format(ac.settings['run']['polygon'])
                    ac.settings['run']['polygon'] = '{}'.format(polygon_new)
                    ac.settings['run']['polygon_clip'] = True
                except:
                    if ac.settings['run']['verbosity'] > 1: print('Provided polygon is not a valid WKT polygon')
                    if ac.settings['run']['verbosity'] > 1: print(ac.settings['run']['polygon'])
                    ac.settings['run']['polygon'] = None
                    pass

            ## read the polygon file
            if os.path.exists(ac.settings['run']['polygon']) & ac.settings['run']['polygon_limit']:
                try:
                    limit = ac.shared.polygon_limit(ac.settings['run']['polygon'])
                    ac.settings['run']['limit'] = [l for l in limit]
                    ac.settings['run']['polygon_clip'] = True
                    if ac.settings['run']['verbosity'] > 1: print('Using limit from polygon envelope: {}'.format(', '.join(['{:}'.format(v) for v in limit])))
                except:
                    if ac.settings['run']['verbosity'] > 1: print('Failed to import polygon {}'.format(ac.settings['run']['polygon']))
        else:
            ac.settings['run']['polygon_clip'] = False

    ## create limit based on station_lon, station_lat, station_box
    if (ac.settings['run']['station_lon'] is not None) &\
       (ac.settings['run']['station_lat'] is not None) &\
       (ac.settings['run']['station_box_size'] is not None) &\
       (ac.settings['run']['limit'] is None):
       site_lat = ac.settings['run']['station_lat']
       site_lon = ac.settings['run']['station_lon']
       box_size = ac.settings['run']['station_box_size']
       print('Creating new limit for position {}N, {}E, box size {} {}'.format(site_lat, site_lon, box_size, ac.settings['run']['station_box_units']))
       if ac.settings['run']['station_box_units'][0] in ['k', 'm']:
           if ac.settings['run']['station_box_units'][0] == 'm': box_size /= 1000.
           ## get approximate distance per degree lon/lat
           dlon, dlat = ac.shared.distance_in_ll(site_lat)
           if type(box_size) is list:
               lat_off, lon_off = (box_size[0]/dlat)/2, (box_size[1]/dlon)/2
           else:
               lat_off, lon_off = (box_size/dlat)/2, (box_size/dlon)/2
       elif ac.settings['run']['station_box_units'][0] in ['d']:
            if type(box_size) is list:
                lat_off, lon_off = box_size[0]/2, box_size[1]/2
            else:
                lat_off, lon_off = box_size/2, box_size/2
       else:
            print('station_box_units={} not configured'.format(ac.settings['run']['station_box_units']))
            return
       ## set new limit
       ac.settings['run']['limit'] = [site_lat-lat_off, site_lon-lon_off, site_lat+lat_off, site_lon+lon_off]
       print('New limit: {}'.format(', '.join(['{:}'.format(v) for v in ac.settings['run']['limit']])))
    ## end create limit based on station information

    ## check limit and add buffer if needed
    if 'limit' in ac.settings['run']:
        if ac.settings['run']['limit'] is not None:
            if len(ac.settings['run']['limit']) != 4:
                print('ROI limit should be four elements in decimal degrees: limit=S,W,N,E')
                print('Provided in the settings:', ac.settings['run']['limit'])
                return

        ## add limit buffer
        if ('limit_buffer' in ac.settings['run']) & (ac.settings['run']['limit'] is not None):
            if (ac.settings['run']['limit_buffer'] is not None):
                if ac.settings['run']['verbosity'] > 1: print('Applying limit buffer of {} {}'.format(ac.settings['run']['limit_buffer'], ac.settings['run']['limit_buffer_units']))
                ac.settings['run']['limit_old'] = [l for l in ac.settings['run']['limit']]

                if ac.settings['run']['limit_buffer_units'][0].lower() == 'd':
                    limit_factor = 1.0, 1.0
                elif ac.settings['run']['limit_buffer_units'][0].lower() in ['m','k']:
                    mean_lat = (ac.settings['run']['limit'][0] + ac.settings['run']['limit'][2]) / 2.
                    dlon, dlat = ac.shared.distance_in_ll(lat=mean_lat)
                    limit_factor = 1/dlat, 1/dlon
                    if ac.settings['run']['limit_buffer_units'][0].lower() == 'm':
                        limit_factor = limit_factor[0]/1000, limit_factor[1]/1000
                else:
                    print('limit_buffer_units={} not configured'.format(ac.settings['run']['limit_buffer_units']))
                    return

                ## compute limit buffer
                limit_buffer = ac.settings['run']['limit_buffer'] * limit_factor[0], \
                               ac.settings['run']['limit_buffer'] * limit_factor[1]
                ac.settings['run']['limit'] = [ac.settings['run']['limit_old'][0] - limit_buffer[0], \
                                               ac.settings['run']['limit_old'][1] - limit_buffer[1], \
                                               ac.settings['run']['limit_old'][2] + limit_buffer[0], \
                                               ac.settings['run']['limit_old'][3] + limit_buffer[1]]
                if ac.settings['run']['verbosity'] > 1:
                    print('Old limit: {}'.format(', '.join(['{:}'.format(v) for v in ac.settings['run']['limit_old']])))
                    print('New limit: {}'.format(', '.join(['{:}'.format(v) for v in ac.settings['run']['limit']])))
        ## set new limits to user settings
        for k in ['limit', 'limit_old']:
            if k in ac.settings['run']: ac.settings['user'][k] = ac.settings['run'][k]
    ## end checking limit settings

    ## make list of lists to process, one list if merging tiles
    inputfile_list = ac.acolite.inputfile_test(ac.settings['run']['inputfile'])

    ## check if tiles need to be merged
    if 'merge_tiles' in ac.settings['run']:
        if ac.settings['run']['merge_tiles']:
            ## figure out whether multiple sets of tiles are given for merging
            ## e.g. through a text inputfile with multiple lines of comma separated scenes
            inputfile = [[]]
            for i in inputfile_list:
                if type(i) == list:
                    inputfile.append(i)
                else:
                    inputfile[0].append(i)
            inputfile_list = [i for i in inputfile if len(i) > 0]
    nruns = len(inputfile_list)

    ## which rgb maps to output from l1r, l2r, l2w files
    l1r_rgb = {'rgb_{}'.format(k): ac.settings['run']['rgb_{}'.format(k)] \
                for k in ac.settings['run']['l1r_rgb_keys'] if 'rgb_{}'.format(k) in ac.settings['run']}
    l1r_rgb_create = any([l1r_rgb[k] for k in l1r_rgb])
    l2r_rgb = {'rgb_{}'.format(k): ac.settings['run']['rgb_{}'.format(k)] \
                for k in ac.settings['run']['l2r_rgb_keys'] if 'rgb_{}'.format(k) in ac.settings['run']}
    l2r_rgb_create = any([l2r_rgb[k] for k in l2r_rgb])
    l2w_rgb = {'rgb_{}'.format(k): ac.settings['run']['rgb_{}'.format(k)] \
                for k in ac.settings['run']['l2w_rgb_keys'] if 'rgb_{}'.format(k) in ac.settings['run']}
    l2w_rgb_create = any([l2w_rgb[k] for k in l2w_rgb])

    ## track processed scenes
    processed = {}
    ## run through bundles to process
    for ni in range(nruns):
        ## bundle to process
        bundle = inputfile_list[ni]
        processed[ni] = {'input': bundle}

        ## save user settings
        settings_file = '{}/acolite_run_{}_l1r_settings_user.txt'.format(ac.settings['run']['output'], ac.settings['run']['runid'])
        ac.acolite.settings.write(settings_file, ac.settings['user'])

        ## run l1 convert
        ret = ac.acolite.acolite_l1r(bundle)
        if ret is None: continue
        if len(ret) == 0: continue
        if len(ret[0]) == 0: continue

        l1r_files, _, l1_bundle = ret
        processed[ni]['l1r'] = l1r_files

        ## project before a/c
        try:
            project = ac.settings['run']['output_projection']
        except:
            project = False
        if (project) & (ac.settings['run']['reproject_before_ac']):
            rep = []
            for ncf in processed[ni]['l1r']:
                ncfo = ac.output.project_acolite_netcdf(ncf)
                if ncfo is None: continue
                rep.append(ncfo)
            if len(rep) > 0:
                processed[ni]['l1r_swath'] = [ncf for ncf in processed[ni]['l1r']]
                l1r_files = [ncf for ncf in rep]
                processed[ni]['l1r'] = l1r_files

        ## save all used settings
        settings_file = '{}/acolite_run_{}_l1r_settings.txt'.format(ac.settings['run']['output'], ac.settings['run']['runid'])
        ac.acolite.settings.write(settings_file, ac.settings['run'])

        ## do atmospheric correction
        l2r_files, l2t_files = [], []
        l2w_files = []
        for l1r in l1r_files:
            gatts = ac.shared.nc_gatts(l1r)
            if 'acolite_file_type' not in gatts: gatts['acolite_file_type'] = 'L1R'
            if ac.settings['run']['l1r_crop']:
                l1r_cropped = ac.output.crop_acolite_netcdf(l1r, output = ac.settings['run']['output'], limit = ac.settings['run']['limit'])
                if l1r_cropped is not None:
                    l1r = l1r_cropped
                else:
                    print('Cropping L1R {} not successful'.format(l1r))
                    continue

            ## L1R geotiff outputs
            if ac.settings['run']['l1r_export_geotiff']: ac.output.nc_to_geotiff(l1r)
            if ac.settings['run']['l1r_export_geotiff_rgb']: ac.output.nc_to_geotiff_rgb(l1r, rgb_datasets = [k.replace('rgb_','') for k in l1r_rgb if l1r_rgb[k]])

            ## rhot RGB
            if (ac.settings['run']['map_l1r']) | (l1r_rgb_create):
                ac.acolite.acolite_map(l1r, plot_all = ac.settings['run']['map_l1r'], settings = l1r_rgb)

            ## do VIS-SWIR atmospheric correction
            if ac.settings['run']['atmospheric_correction']:
                if gatts['acolite_file_type'].startswith('L1R'):
                    ## run dsf or exp
                    if (ac.settings['run']['atmospheric_correction_method'] in ['dark_spectrum', 'exponential']):
                        ret = ac.acolite.acolite_l2r(l1r)
                        if ret is None:
                            l2r = []
                        elif len(ret) != 2:
                            l2r = []
                        else:
                            l2r, _ = ret
                    ## run radcor
                    elif (ac.settings['run']['atmospheric_correction_method'] == 'radcor'):
                        l2r = ac.adjacency.radcor.radcor(l1r)
                        if l2r is None: l2r = []
                    else:
                        print('Option atmospheric_correction_method={} not configured, use dark_spectrum, radcor, or exponential'.format(ac.settings['run']['atmospheric_correction_method']))
                else:
                    l2r = '{}'.format(l1r)

                ## run glad after dsf (not recommended as a general method)
                if (ac.settings['run']['adjacency_correction']) &\
                    (ac.settings['run']['adjacency_correction_method'] == 'glad') &\
                    (ac.settings['run']['atmospheric_correction_method'] == 'dark_spectrum'):
                    if len(l2r) > 0:
                        ret = ac.adjacency.glad.glad_l2r(l2r, settings = ac.settings['run'], verbosity = ac.config['verbosity'])
                        l2r = [] if ret is None else ret

                ## if we have multiple l2r files
                if (len(l2r) > 0):
                    if type(l2r) is not list: l2r = [l2r]
                    l2r_files+=l2r

                    ## run pansharpening
                    if ac.settings['run']['pans']:
                        for ncf in l2r:
                            pr = ac.acolite.acolite_pans(ncf)
                            if pr is not None:
                                if 'l2r_pans' not in processed[ni]: processed[ni]['l2r_pans']=[]
                                processed[ni]['l2r_pans'].append(pr)

                    ## run outputs
                    l2r_ = [ncf for ncf in l2r]
                    if 'l2r_pans' in processed[ni]: l2r_ += processed[ni]['l2r_pans'] ## include pansharpened L2R
                    for ncf in l2r_:
                        ## L2R geotiff outputs
                        if ac.settings['run']['l2r_export_geotiff']: ac.output.nc_to_geotiff(ncf)
                        if ac.settings['run']['l2r_export_geotiff_rgb']: ac.output.nc_to_geotiff_rgb(ncf, rgb_datasets = [k.replace('rgb_','') for k in l2r_rgb if l2r_rgb[k]])

                        ## make rgb rhos maps
                        if (ac.settings['run']['map_l2r']) | l2r_rgb_create:
                            ac.acolite.acolite_map(ncf, plot_all = ac.settings['run']['map_l2r'], settings = l2r_rgb)

                    ## compute l2w parameters
                    if ac.settings['run']['l2w_parameters'] is not None:
                        for ncf in l2r:
                            ret = ac.acolite.acolite_l2w(ncf)
                            if ret is not None:
                                if ac.settings['run']['l2w_export_geotiff']: ac.output.nc_to_geotiff(ret)
                                l2w_files.append(ret)

                                ## make l2w maps
                                if (ac.settings['run']['map_l2w']) | l2w_rgb_create:
                                    ac.acolite.acolite_map(ret, plot_all=ac.settings['run']['map_l2w'], settings =l2w_rgb)

            ## run TACT thermal atmospheric correction
            if ac.settings['run']['tact_run']:
                ret = ac.tact.tact_gem(l1r, settings = ac.settings['run'])
                if ret is not None:
                    l2t_files.append(ret)
                    if ac.settings['run']['l2t_export_geotiff']: ac.output.nc_to_geotiff(ret)

                    ## make l2t maps
                    if ac.settings['run']['map_l2t']: ac.acolite.acolite_map(ret, plot_all=ac.settings['run']['map_l2t'])

        if len(l2r_files) > 0: processed[ni]['l2r'] = l2r_files
        if len(l2t_files) > 0: processed[ni]['l2t'] = l2t_files
        if len(l2w_files) > 0: processed[ni]['l2w'] = l2w_files

    ## reproject data
    try:
        project = (ac.settings['run']['output_projection']) & (~ac.settings['run']['reproject_before_ac'])
    except:
        project = False
    if project:
        for output in ac.settings['run']['reproject_outputs']:
            otype = output.lower()
            pkeys = processed.keys()
            for i in pkeys:
                reprojected = []
                if otype not in processed[i]: continue
                for ncf in processed[i][otype]:
                    ncfo = ac.output.project_acolite_netcdf(ncf)
                    if ncfo is None: continue
                    reprojected.append(ncfo)

                    ## make rgb  maps
                    if (otype == 'l1r') & (l1r_rgb_create | ac.settings['run']['map_l1r']):
                        ac.acolite.acolite_map(ncfo, plot_all = ac.settings['run']['map_l1r'], settings = l1r_rgb)
                        rgb_datasets = [k.replace('rgb_','') for k in l1r_rgb if l1r_rgb[k]]

                    ## make rgb  maps
                    if (otype == 'l2r') & (l2r_rgb_create | ac.settings['run']['map_l2r']):
                        ac.acolite.acolite_map(ncfo, plot_all = ac.settings['run']['map_l2r'], settings = l2r_rgb)
                        rgb_datasets = [k.replace('rgb_','') for k in l2r_rgb if l2r_rgb[k]]

                    ## make rgb and other maps
                    if (otype == 'l2w') & (l2w_rgb_create | ac.settings['run']['map_l2w']):
                        ac.acolite.acolite_map(ncfo, plot_all = ac.settings['run']['map_l2w'], settings = l2w_rgb)
                        rgb_datasets = [k.replace('rgb_','') for k in l2w_rgb if l2w_rgb[k]]

                    ## make tact maps
                    if (otype == 'l2t') & (ac.settings['run']['map_l2t']):
                        ac.acolite.acolite_map(ncfo, plot_all = ac.settings['run']['map_l2t'])

                    ## output geotiffs
                    if '{}_export_geotiff'.format(otype) in ac.settings['run']:
                        if ac.settings['run']['{}_export_geotiff'.format(otype)]:
                            ac.output.nc_to_geotiff(ncfo)

                    ## output rgb geotiff
                    if '{}_export_geotiff_rgb'.format(otype) in ac.settings['run']:
                        if ac.settings['run']['{}_export_geotiff_rgb'.format(otype)]:
                            ac.output.nc_to_geotiff_rgb(ncfo, rgb_datasets = rgb_datasets)
                processed[i]['{}_reprojected'.format(otype)] = reprojected
    ## end reproject data

    ## end processing loop
    log.__del__()

    ## remove files
    for ni in processed:
        ## remove output netcdfs
        for level in ['l1r', 'l2r', 'l2r_pans', 'l2t', 'l2w']:
            if level not in processed[ni]: continue
            if '{}_delete_netcdf'.format(level) not in ac.settings['run']: continue
            if ac.settings['run']['{}_delete_netcdf'.format(level)]:
                ## run through images and delete them
                for f in processed[ni][level]:
                    if os.path.exists(f):
                        print('Deleting {}'.format(f))
                        os.remove(f)
                    ## also delete pan file if it exists
                    if level == 'l1r':
                        panf = f.replace('_L1R.nc', '_L1R_pan.nc')
                        if os.path.exists(panf):
                            print('Deleting {}'.format(panf))
                            os.remove(panf)

    ## remove log and settings files for this run
    try:
        if ac.settings['run']['delete_acolite_run_text_files']:
            tfiles = glob.glob('{}/acolite_run_{}_*.txt'.format(ac.settings['run']['output'], ac.settings['run']['runid']))
            for tf in tfiles: os.remove(tf)
            if 'polygon_old' in ac.settings['run']:
                if ac.settings['run']['polygon_old'] != ac.settings['run']['polygon']:
                    os.remove(ac.settings['run']['polygon'])
    except KeyError:
        print('Not removing text files as "delete_acolite_run_text_files" is not in settings.')
    except BaseException as err:
        print('Could not remove ACOLITE text files in directory {}'.format(ac.settings['run']['output']))

    ## remove output directory if delete_acolite_output_directory is True and the directory is empty
    try:
        if (ac.settings['run']['delete_acolite_output_directory']) & (new_path is not None):
            output = os.path.abspath(ac.settings['run']['output'])
            output_split = output.split(os.path.sep)
            ## create list of directory levels to delete
            to_delete = []
            test_path = ''
            for l in output_split:
                test_path += l+os.path.sep
                if len(test_path) < len(new_path): continue
                to_delete.append(test_path)
            ## start from last directory level
            to_delete.reverse()
            for td in to_delete:
                os.rmdir(td)
                print('Removed {}'.format(td))
    except KeyError:
        print('Not testing output directory as "delete_acolite_output_directory" is not in settings.')
    except OSError as err:
        print('Could not remove output directory {}'.format(ac.settings['run']['output']))
        if (err.errno == 66): print('Output directory not empty {}'.format(ac.settings['run']['output']))

    return(processed)
