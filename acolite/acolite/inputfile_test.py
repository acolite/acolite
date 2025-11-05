## def inputfile_test
## function to test/parse inputfile given to acolite
## written by Quinten Vanhellemont, RBINS
## 2022-03-04
## modifications: 2022-09-05 (QV) remove trailing slash from file paths in txt file
##                2022-11-14 (QV) changed appending files from txt file to inputfile list
##                2023-09-12 (QV) added attempt to download from CDSE
##                2023-09-21 (QV) added attempt to download from EarthExplorer
##                2023-10-21 (QV) added ECOSTRESS download from EarthExplorer
##                2024-04-27 (QV) moved APIs
##                2025-01-23 (QV) test if scene exists in download directory, added S2C
##                2025-02-05 (QV) add oceandata directdataaccess for MERIS
##                2025-02-09 (QV) added config download directory
##                2025-05-13 (QV) test if inputfile exists, then spit on semicolon, then on comma

def inputfile_test(inputfile):
    import os, mimetypes
    import acolite as ac

    ## check if a list of files is given
    if type(inputfile) == str:
        if os.path.exists(inputfile):
            tmp_files = [inputfile]
        elif inputfile.startswith('https:') & inputfile.endswith('.zarr'):
            tmp_files = [inputfile]
        elif ';' in inputfile:
            tmp_files = inputfile.split(';')
        else:
            tmp_files = inputfile.split(',')
    elif type(inputfile) == list:
        tmp_files = inputfile
    else:
        if ac.config['verbosity'] > 0: print('Inputfile {} not recognised'.format(inputfile))
        tmp_files = []

    ## run through files
    inputfile_list = []
    for file in tmp_files:
        if len(file) == 0: continue
        file = file.strip() ## strip spaces

        if file.startswith('https:') & file.endswith('.zarr'):
            if ac.config['verbosity'] > 0: print('Assuming {} is zarr url.'.format(file))
            inputfile_list.append(file)
            continue
            
        if not os.path.exists(file):
            if ac.config['verbosity'] > 0: print('Path {} does not exist.'.format(file))
            ## try and download from CDSE or EarthExplorer
            if ac.settings['run']['scene_download']:
                if ac.settings['run']['scene_download_directory'] is not None:
                    ddir = ac.settings['run']['scene_download_directory']
                elif ac.config['scene_download_directory'] is not None:
                    ddir = ac.config['scene_download_directory']
                else:
                    ddir = ac.settings['run']['output']
                ## find out data source to use
                bn = os.path.basename(file)
                local_file = '{}/{}'.format(ddir, bn)
                ## test if file exists in download dir
                if os.path.exists(local_file):
                    file = '{}'.format(local_file)
                    print('Scene exists at {}'.format(file))
                ## otherwise try and download it
                else:
                    if bn[0:3] in ['S2A', 'S2B', 'S2C', 'S3A', 'S3B']:
                        download_source = 'CDSE'
                    elif bn[0:4] in ['LC08', 'LO08', 'LT08', 'LC09', 'LO09', 'LT09', 'LT04', 'LT05', 'LE07']:
                        download_source = 'EarthExplorer'
                    elif 'ECOSTRESS' in bn:
                        download_source = 'EarthExplorer'
                    elif bn.startswith('EN1_MDSI_MER_'): ## MERIS data from oceandata directdataaccess
                        download_source = 'oceandata'
                        if not bn.endswith('.ZIP'): bn += '.ZIP'
                    elif bn.startswith('H2') & ('.L1B_ISS' in bn): ## HICO
                        download_source = 'oceandata'
                        #if not bn.endswith('.bz2'): bn += '.bz2'
                    else:
                        print('Could not identify download source for scene {}'.format(file))
                        continue

                    ## update local file if bn was updated
                    local_file = '{}/{}'.format(ddir, bn)

                    if ac.config['verbosity'] > 0: print('Attempting download of scene {} from {}.'.format(file, download_source))
                    if ac.config['verbosity'] > 0: print('Target directory {}'.format(ddir))
                    if ac.config['verbosity'] > 0: print('Querying {}'.format(download_source))

                    ## Copernicus Data Space Ecosystem
                    if download_source == 'CDSE':
                        urls, scenes = ac.api.cdse.query(scene=bn)
                        if ac.config['verbosity'] > 0: print('Downloading from {}'.format(download_source))
                        local_scenes = ac.api.cdse.download(urls, output = ddir, verbosity = ac.config['verbosity'])
                    ## EarthExplorer
                    if download_source == 'EarthExplorer':
                        entity_list, identifier_list, dataset_list = ac.api.earthexplorer.query(scene=bn)
                        if ac.config['verbosity'] > 0: print('Downloading from {}'.format(download_source))
                        local_scenes = ac.api.earthexplorer.download(entity_list, dataset_list, identifier_list, output = ddir, verbosity = ac.config['verbosity'])
                    ## oceandata directdataaccess
                    if download_source == 'oceandata':
                        local_scenes = []
                        url  = '{}/{}'.format(ac.config['oceandata_url'], bn)
                        if ac.config['verbosity'] > 0: print('Downloading from {}'.format(download_source))
                        try:
                            ac.shared.download_file(url, local_file)
                        except:
                            if ac.config['verbosity'] > 0: print('Could not complete download from {}'.format(download_source))
                        if os.path.exists(local_file): local_scenes.append(local_file)
                    if len(local_scenes) == 1:  file = '{}'.format(local_scenes[0])
                    if not os.path.exists(file): continue
            else:
                continue

        ##  remove trailing slash
        if file[-1] == os.sep: file = file[0:-1]

        if os.path.isdir(file):
            inputfile_list.append(file)
        else:
            mime = mimetypes.guess_type(file)
            if mime[0] != 'text/plain':
                if os.path.exists(file): inputfile_list.append(file) ## assume we can process this file
                continue
            with open(file, 'r') as f:
                for l in f.readlines():
                    l = l.strip()
                    if len(l) == 0: continue
                    cfiles = []
                    for fn in l.split(','):
                        fn = fn.strip()
                        if fn[-1] == os.sep: fn = fn[0:-1]
                        if os.path.exists(fn):
                            cfiles.append(fn)
                        else:
                            if ac.config['verbosity'] > 0: print('Path {} does not exist.'.format(fn))
                    if len(cfiles)>0: inputfile_list += cfiles
    return(inputfile_list)
