## def inputfile_test
## function to test/parse inputfile given to acolite
## written by Quinten Vanhellemont, RBINS
## 2022-03-04
## modifications:

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
        if not os.path.exists(file):
            if ac.config['verbosity'] > 0: print('Path {} does not exist.'.format(file))
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
                        if os.path.exists(fn):
                            cfiles.append(fn)
                        else:
                            if ac.config['verbosity'] > 0: print('Path {} does not exist.'.format(fn))
                    if len(cfiles)>0: inputfile_list.append(cfiles)
    return(inputfile_list)
