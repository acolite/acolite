## def acolite_run
## runs acolite processing for given settings file/dict, new version
## written by Quinten Vanhellemont, RBINS
## 2021-04-01
## modifications: 2021-04-14 (QV) added output to settings if not configured
##                2021-04-15 (QV) test/parse input files
##                2022-03-04 (QV) moved inputfile testing to inputfile_test

def acolite_run(settings, inputfile=None, output=None):
    import glob, datetime, os
    import acolite as ac

    print('Running ACOLITE processing - {}'.format(ac.version))
    ## time of processing start
    time_start = datetime.datetime.now()

    ## get user settings
    ## these are updated with sensor specific settings in acolite_l2r
    setu = ac.acolite.settings.parse(None, settings=settings, merge=False)
    if 'runid' not in setu: setu['runid'] = time_start.strftime('%Y%m%d_%H%M%S')
    if 'output' not in setu: setu['output'] = os.getcwd()
    if 'verbosity' in setu: ac.config['verbosity'] = int(setu['verbosity'])

    ## workaround for outputting rhorc and bt
    if 'l2w_parameters' in setu:
        if setu['l2w_parameters'] is not None:
            for par in setu['l2w_parameters']:
                if 'rhorc' in par: setu['output_rhorc'] = True
                if 'bt' == par[0:2]: setu['output_bt'] = True

    ## log file for l1r generation
    log_file = '{}/acolite_run_{}_log_file.txt'.format(setu['output'],setu['runid'])
    log = ac.acolite.logging.LogTee(log_file)
    print('Run ID - {}'.format(setu['runid']))

    ## earthdata credentials from settings file
    for k in ['EARTHDATA_u', 'EARTHDATA_p']:
        kv = setu[k] if k in setu else ac.config[k]
        if len(kv) == 0: continue
        os.environ[k] = kv

    ## set settings from launch_acolite
    if inputfile is not None: setu['inputfile'] = inputfile
    if output is not None: setu['output'] = output

    ## check if we have anything to do
    if 'inputfile' not in setu:
        print('Nothing to do. Did you provide a settings file or inputfile? Exiting.')
        return()
    else:
        nscenes = len(setu['inputfile'])

    ## get defaults settings for l1r processing
    setu_l1r = ac.acolite.settings.parse(None, settings=setu, merge=False)

    ## make list of lists to process, one list if merging tiles
    inputfile_list = ac.acolite.inputfile_test(setu['inputfile'])

    ## check if tiles need to be merged
    if 'merge_tiles' in setu:
        if setu['merge_tiles']:
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

    ## track processed scenes
    processed = {}
    ## run through bundles to process
    for ni in range(nruns):
        ## bundle to process
        bundle = inputfile_list[ni]
        processed[ni] = {'input': bundle}

        ## save user settings
        settings_file = '{}/acolite_run_{}_l1r_settings_user.txt'.format(setu_l1r['output'],setu_l1r['runid'])
        ac.acolite.settings.write(settings_file, setu_l1r)

        ## run l1 convert
        ret = ac.acolite.acolite_l1r(bundle, setu_l1r)
        if len(ret) == 0: continue
        l1r_files, l1r_setu = ret
        processed[ni]['l1r'] = l1r_files

        ## save all used settings
        settings_file = '{}/acolite_run_{}_l1r_settings.txt'.format(l1r_setu['output'],l1r_setu['runid'])
        ac.acolite.settings.write(settings_file, l1r_setu)

        ## do atmospheric correction
        l2r_files, l2t_files = [], []
        l2w_files = []
        for l1r in l1r_files:
            gatts = ac.shared.nc_gatts(l1r)
            if 'acolite_file_type' not in gatts: gatts['acolite_file_type'] = 'L1R'
            if l1r_setu['l1r_export_geotiff']: ac.output.nc_to_geotiff(l1r, match_file = l1r_setu['export_geotiff_match_file'],
                                                            cloud_optimized_geotiff = l1r_setu['export_cloud_optimized_geotiff'],
                                                            skip_geo = l1r_setu['export_geotiff_coordinates'] is False)
            if l1r_setu['l1r_export_geotiff_rgb']: ac.output.nc_to_geotiff_rgb(l1r, settings = l1r_setu)

            ## rhot RGB
            if l1r_setu['rgb_rhot']:
                l1r_setu_ = {k: l1r_setu[k] for k in l1r_setu}
                l1r_setu_['rgb_rhos'] = False
                ac.acolite.acolite_map(l1r, settings = l1r_setu_, plot_all=False)

            ## do VIS-SWIR atmospheric correction
            if l1r_setu['atmospheric_correction']:
                if gatts['acolite_file_type'] == 'L1R':
                    ## run ACOLITE
                    ret = ac.acolite.acolite_l2r(l1r, settings = l1r_setu, verbosity = ac.config['verbosity'])
                    if len(ret) != 2: continue
                    l2r, l2r_setu = ret
                else:
                    l2r = '{}'.format(l1r)
                    l2r_setu = ac.acolite.settings.parse(gatts['sensor'], settings=l1r_setu)

                if (l2r_setu['adjacency_correction']):
                    ret = None
                    ## acstar3 adjacency correction
                    if (l2r_setu['adjacency_method']=='acstar3'):
                        ret = ac.adjacency.acstar3.acstar3(l2r, setu = l2r_setu, verbosity = ac.config['verbosity'])
                    ## GLAD
                    if (l2r_setu['adjacency_method']=='glad'):
                        ret = ac.adjacency.glad.glad_l2r(l2r, verbosity = ac.config['verbosity'], settings=l2r_setu)
                    l2r = [] if ret is None else ret

                ## if we have multiple l2r files
                if type(l2r) is not list: l2r = [l2r]
                l2r_files+=l2r

                if l2r_setu['l2r_export_geotiff']:
                    for ncf in l2r:
                        ac.output.nc_to_geotiff(ncf, match_file = l2r_setu['export_geotiff_match_file'],
                                                cloud_optimized_geotiff = l1r_setu['export_cloud_optimized_geotiff'],
                                                skip_geo = l2r_setu['export_geotiff_coordinates'] is False)

                        if l2r_setu['l2r_export_geotiff_rgb']: ac.output.nc_to_geotiff_rgb(ncf, settings = l2r_setu)


                ## make rgb rhos maps
                if l2r_setu['rgb_rhos']:
                    l2r_setu_ = {k: l1r_setu[k] for k in l2r_setu}
                    l2r_setu_['rgb_rhot'] = False
                    for ncf in l2r:
                        ac.acolite.acolite_map(ncf, settings = l2r_setu_, plot_all=False)

                ## compute l2w parameters
                if l2r_setu['l2w_parameters'] is not None:
                    if type(l2r_setu['l2w_parameters']) is not list: l2r_setu['l2w_parameters'] = [l2r_setu['l2w_parameters']]
                    for ncf in l2r:
                        ret = ac.acolite.acolite_l2w(ncf, settings=l2r_setu)
                        if ret is not None:
                            if l2r_setu['l2w_export_geotiff']: ac.output.nc_to_geotiff(ret, match_file = l2r_setu['export_geotiff_match_file'],
                                                                            cloud_optimized_geotiff = l1r_setu['export_cloud_optimized_geotiff'],
                                                                            skip_geo = l2r_setu['export_geotiff_coordinates'] is False)
                            l2w_files.append(ret)

                            ## make l2w maps
                            if l2r_setu['map_l2w']:
                                ac.acolite.acolite_map(ret, settings=l2r_setu)
                            ## make l2w rgb
                            if l2r_setu['rgb_rhow']:
                                l2r_setu_ = {k: l1r_setu[k] for k in l2r_setu}
                                l2r_setu_['rgb_rhot'] = False
                                l2r_setu_['rgb_rhos'] = False
                                ac.acolite.acolite_map(ret, settings=l2r_setu_, plot_all=False)

            ## run TACT thermal atmospheric correction
            if l1r_setu['tact_run']:
                ret = ac.tact.tact_gem(l1r, output = l1r_setu['output'], verbosity=ac.config['verbosity'],
                                            output_atmosphere = l1r_setu['tact_output_atmosphere'],
                                            output_intermediate = l1r_setu['tact_output_intermediate'])
                if ret != ():
                    l2t_files.append(ret)
                    if l1r_setu['l2t_export_geotiff']: ac.output.nc_to_geotiff(ret, match_file = l1r_setu['export_geotiff_match_file'],
                                                                               cloud_optimized_geotiff = l1r_setu['export_cloud_optimized_geotiff'],
                                                                               skip_geo = l1r_setu['export_geotiff_coordinates'] is False)

                    ## make l2t maps
                    if l1r_setu['tact_map']: ac.acolite.acolite_map(ret, settings=l1r_setu)

        if len(l2r_files) > 0: processed[ni]['l2r'] = l2r_files
        if len(l2t_files) > 0: processed[ni]['l2t'] = l2t_files
        if len(l2w_files) > 0: processed[ni]['l2w'] = l2w_files

    ## reproject data
    try:
        project = l1r_setu['output_projection']
    except:
        project = False
    if project:
        for output in l1r_setu['reproject_outputs']:
            otype = output.lower()
            pkeys = processed.keys()
            for i in pkeys:
                reprojected = []
                if otype not in processed[i]: continue
                for ncf in processed[i][otype]:
                    ncfo = ac.output.project_acolite_netcdf(ncf, settings=settings)
                    if ncfo == (): continue
                    reprojected.append(ncfo)

                    ## make rgb  maps
                    if (otype == 'l1r') & (l1r_setu['rgb_rhot']):
                        l1r_setu_ = {k: l1r_setu[k] for k in l1r_setu}
                        l1r_setu_['rgb_rhos'] = False
                        ac.acolite.acolite_map(ncfo, settings = l1r_setu_, plot_all=False)
                    ## make rgb  maps
                    if (otype == 'l2r') & (l1r_setu['rgb_rhos']):
                        l1r_setu_ = {k: l1r_setu[k] for k in l1r_setu}
                        l1r_setu_['rgb_rhot'] = False
                        ac.acolite.acolite_map(ncfo, settings = l1r_setu_, plot_all=False)

                    ## make other  maps
                    if (otype == 'l2w') & (l1r_setu['map_l2w']):
                        l1r_setu_ = {k: l1r_setu[k] for k in l1r_setu}
                        l1r_setu_['rgb_rhos'] = False
                        l1r_setu_['rgb_rhot'] = False
                        ac.acolite.acolite_map(ncfo, settings = l1r_setu_, plot_all=True)

                    ## output geotiffs
                    if '{}_export_geotiff'.format(otype) in l1r_setu:
                        if l1r_setu['{}_export_geotiff'.format(otype)]:
                            ac.output.nc_to_geotiff(ncfo, match_file = l1r_setu['export_geotiff_match_file'],
                                                                        cloud_optimized_geotiff = l1r_setu['export_cloud_optimized_geotiff'],
                                                                        skip_geo = l1r_setu['export_geotiff_coordinates'] is False)
                    ## output rgb geotiff
                    if '{}_export_geotiff_rgb'.format(otype) in l1r_setu:
                        if l1r_setu['{}_export_geotiff_rgb'.format(otype)]:
                            ac.output.nc_to_geotiff_rgb(ncfo, settings = l1r_setu)
                processed[i]['{}_reprojected'.format(otype)] = reprojected
    ## end reproject data

    ## end processing loop
    log.__del__()

    ## remove files
    for ni in processed:
        for level in ['l1r', 'l2r', 'l2t', 'l2w']:
            if level not in processed[ni]: continue
            if '{}_delete_netcdf'.format(level) not in l1r_setu: continue
            if l1r_setu['{}_delete_netcdf'.format(level)]:
                ## run through images and delete them
                for f in processed[ni][level]:
                    os.remove(f)
                    ## also delete pan file if it exists
                    if level == 'l1r':
                        panf = f.replace('_L1R.nc', '_L1R_pan.nc')
                        if os.path.exists(panf):
                            os.remove(panf)
                ## replace by empty list
                processed[ni][level] = []

    ## remove log and settings files
    try:
        delete_text = l1r_setu['delete_acolite_run_text_files']
        op, ri = l1r_setu['output'], l1r_setu['runid']
    except:
        delete_text = False if 'delete_acolite_run_text_files' not in setu_l1r else setu_l1r['delete_acolite_run_text_files']
        op, ri = setu_l1r['output'], setu_l1r['runid']

    if delete_text:
        tfiles = glob.glob('{}/acolite_run_{}_*.txt'.format(op, ri))
        for tf in tfiles:
            os.remove(tf)

    return(processed)
