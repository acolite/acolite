## def inputfile_test
## function to test/parse inputfile given to acolite
## written by Quinten Vanhellemont, RBINS
## 2022-03-04
## modifications: 2022-09-05 (QV) remove trailing slash from file paths in txt file
##                2022-11-14 (QV) changed appending files from txt file to inputfile list
##                2023-09-12 (QV) added attempt to download from CDSE

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
            ## try and download from CDSE
            if ac.settings['run']['cdse_download']:
                if ac.config['verbosity'] > 0: print('Attempting download of scene {} from CDSE.'.format(file))
                ddir = ac.settings['run']['cdse_download_directory']
                if ddir is None: ddir = ac.settings['run']['output']
                if ac.config['verbosity'] > 0: print('Querying CDSE')
                urls, scenes = ac.cdse.query(scene=os.path.basename(inputfile))
                if ac.config['verbosity'] > 0: print('Downloading from CDSE')
                local_scenes = ac.cdse.download(urls, output = ddir)
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
