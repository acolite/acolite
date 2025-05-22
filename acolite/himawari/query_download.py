## def query_download
## query and download Himawari data using himawari_download_script
## get Full Disk Himawari Standard Data for VSWIR bands B01-B06
## needs your credentials in your netrc file:
##
## machine himawari.diasjp.net
## login YOUR_LOGIN
## password YOUR_PASSWORD
##
## written by Quinten Vanhellemont, RBINS
## 2025-05-22
## modifications:

def query_download(date_start, date_end = None, time_diff = 660,
                   bands = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06'], ## VSWIR bands
                   data_type = 'HS' , area_type = 'FLDK',
                   download_dir = None, himawari_download_script_dir = None, netrc = None,
                   strip_directory_name = False, download = True):

    import dateutil.parser, datetime
    import os, sys, subprocess
    cwd = os.getcwd()

    ## test script dir
    if himawari_download_script_dir is None:
        if 'himawari_download_script_dir' in os.environ:
            himawari_download_script_dir = os.environ['himawari_download_script_dir']
        elif 'himawari_download_script_dir' in ac.config:
            himawari_download_script_dir = ac.config['himawari_download_script_dir']
        else:
            print('Could not determine himawari_download_script_dir.')
            return

    ## download dir
    if (download) & (download_dir is None):
        print('Please provide download_dir')
        return
    if not os.path.exists(download_dir): os.makedirs(download_dir)

    ## test
    if netrc is not None:
        if not os.path.exists(netrc):
            print('Please provide valid netrc file')
            return

    ## set up dates
    if type(date_start) == str:
        date1 = dateutil.parser.parse(date_start)
    elif type(date_start) == datetime.datetime:
        date1 = date_start
    else:
        print('Not a valid start date: {}'.format(date_start))
        return
    if date_end is None:
        date2 = date1 + datetime.timedelta(seconds=time_diff)
    elif type(date_end) == str:
        date2 = dateutil.parser.parse(date_end)
    elif type(date_end) == datetime.datetime:
        date2 = date_end
    else:
        print('Not a valid end date: {}'.format(date_end))
        return

    ## get date name for file testing
    date_name = date1.strftime('%Y%m%d_%H%M')

    ## query command
    cmd = ['{}/himawari-ls.py'.format(himawari_download_script_dir), '-T {}'.format(data_type), '-A {}'.format(area_type),\
          '-f {}'.format(date1.isoformat()[0:16]), '-t {}'.format(date2.isoformat()[0:16])]
    if netrc is not None: cmd += ['-n {}'.format(netrc)]

    ## run query
    p = subprocess.run(' '.join(cmd), shell=True, stdout=subprocess.PIPE)
    if (p.returncode != 0):
        print('Error querying:')
        print(p.stdout.decode(encoding='utf-8').split('\n'))
        return

    ## parse returned files
    files_ = p.stdout.decode(encoding='utf-8').split('\n')
    files_.sort()
    p = None

    ## identify files for requested cycle
    ## query sometimes adds files from previous cycle
    files = []
    for file in files_:
        ## remove newlines
        file = file.strip()
        if len(file) == 0: continue

        bn = os.path.basename(file)
        dn = os.path.dirname(file)

        ## test if file is from the requested cycle
        if (date_end is None) & (bn[7:20] != date_name): continue

        ## test if band in bands
        if bn[21:24] not in bands: continue

        files.append(file)

    if not download: return(files)

    ofiles = []
    for file in files:
        bn = os.path.basename(file)
        dn = os.path.dirname(file)

        if strip_directory_name:
            ofile = '{}/{}'.format(download_dir, bn)
        else:
            ofile = '{}/{}/{}'.format(download_dir, dn, bn)

        if not os.path.exists(ofile):
            if not os.path.exists(os.path.dirname(ofile)):
                os.makedirs(os.path.dirname(ofile))

            print('Downloading {}'.format(file))
            cmd = ['{}/himawari-dl.py'.format(himawari_download_script_dir), file, '-o "{}"'.format(ofile)]
            if netrc is not None: cmd += ['-n {}'.format(netrc)]
            cmd = ' '.join(cmd)
            p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
            if p.returncode == 1:
                print('Error for file {}'.format(file))
                if os.path.exists(ofile): os.remove(ofile)
            else:
                print('Success for file {}'.format(file))
            p = None

        if os.path.exists(ofile): ofiles.append(ofile)
    return(ofiles)        
