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

def inputfile_test(inputfile):
    import os, mimetypes
    import acolite as ac

    ## check if a list of files is given
    if type(inputfile) == str:
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
        if not os.path.exists(file):
            if ac.config['verbosity'] > 0: print('Path {} does not exist.'.format(file))
            ## try and download from CDSE or EarthExplorer
            if ac.settings['run']['scene_download']:
                ddir = ac.settings['run']['scene_download_directory']
                if ddir is None: ddir = ac.settings['run']['output']
                ## find out data source to use
                bn = os.path.basename(inputfile)
                ## test if file exists in download dir
                if os.path.exists('{}/{}'.format(ddir, os.path.basename(file))):
                    file = '{}/{}'.format(ddir, os.path.basename(file))
                    print('Scene exists at {}'.format(file))
                ## otherwise try and download it
                else:
                    if bn[0:3] in ['S2A', 'S2B', 'S2C', 'S3A', 'S3B']:
                        download_source = 'CDSE'
                    elif bn[0:4] in ['LC08', 'LO08', 'LT08', 'LC09', 'LO09', 'LT09', 'LT04', 'LT05', 'LE07']:
                        download_source = 'EarthExplorer'
                    elif 'ECOSTRESS' in bn:
                        download_source = 'EarthExplorer'
                    else:
                        print('Could not identify download source for scene {}'.format(file))
                        continue

                    if ac.config['verbosity'] > 0: print('Attempting download of scene {} from {}.'.format(file, download_source))
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
