## def sftp.upload
## uploads all files in path to sftp, currently only with stat file size check
## written by Quinten Vanhellemont, RBINS
## 2026-04-20
## modifications: 2026-04-20 (QV) added skip paths, log_uploaded_files

def upload(host, local_path, port = 22, remote_path = '/',
            override = False, machine_append = '', skip_paths = [], log_uploaded_files = False):

    import os, glob, paramiko, netrc

    ## get netrc machine definition
    machine = '{}{}'.format(host, machine_append)
    ## get credentials
    try:
        nr = netrc.netrc(os.environ['NETRC'])
    except KeyError:
        nr = netrc.netrc()

    ret = nr.authenticators(machine)
    if ret is not None:
        username, account, password = ret
    else:
        print('Could not load .netrc authentication for machine {}'.format(machine))
        return

    ## check log
    logged_files = []
    if log_uploaded_files:
        if local_path[-1] == os.path.sep:
            log_file = '{}_upload.log'.format(local_path[0:-1])
        else:
            log_file = '{}_upload.log'.format(local_path)
        ## remove log file if override
        if (os.path.exists(log_file)) & (override):
            os.remove(log_file)
            print('Removed log file at {} because override = True'.format(log_file))
        if os.path.exists(log_file):
            ## read log
            with open(log_file, 'r', encoding = 'utf-8') as f:
                for line in f.readlines():
                    line = line.strip()
                    if len(line) > 0:
                        logged_files.append(line)

    ## find local files
    local_files = glob.glob('{}/**'.format(local_path), recursive = True)
    local_files.sort()
    print('Recursively found {} files in {}'.format(len(local_files), local_path))
    if len(local_files) == 0: return

    ## make a list of given skip paths
    if type(skip_paths) != list: skip_paths = [skip_paths]

    ## set up transport
    transport = paramiko.Transport((host, port))
    transport.connect(username = username, password = password)

    ## connect and upload files
    with paramiko.SFTPClient.from_transport(transport) as sftp:
        ## run through files
        for local_file in local_files:
            if os.path.isdir(local_file): continue ## skip directories

            ## get file info and remote path
            dn = os.path.dirname(local_file)
            bn = os.path.basename(local_file)
            local_stat = os.stat(local_file)
            remote_dn = dn.replace(local_path, remote_path + '/')
            remote_file = '{}/{}'.format(remote_dn, bn)

            ## skip file if logged
            if (remote_file in logged_files) & (not override): continue

            ## skip matches in directory name
            skip_file = False
            for skip in skip_paths:
                if type(skip) is not str: continue
                if skip in dn: skip_file = True
            if skip_file: continue

            ## test if file extist and get size
            upload = True
            try:
                remote_stat = sftp.stat(remote_file)
                if (remote_stat.st_size != local_stat.st_size) | (override):
                    sftp.remove(remote_file)
                else:
                    upload = False
            except:
                pass

            ## if we need to upload
            if upload:
                print('Uploading {}'.format(remote_file))
                try:
                    sftp.chdir(remote_dn)
                except:
                    sftp.chdir('/') ## start at base dir
                    remote_path_walk = remote_dn.split('/')
                    # get to target directory
                    for rp in remote_path_walk:
                       try:
                           sftp.chdir(rp)
                       except IOError:
                           print('Creating {}'.format(rp))
                           sftp.mkdir(rp)
                           sftp.chdir(rp)

                sftp.put(local_file, bn)

                ## log upload
                if log_uploaded_files:
                    logged_files.append(remote_file)
                    with open(log_file, 'a', encoding = 'utf-8') as f:
                        f.write(remote_file + '\n')
            else:
                print('File {} exists with correct size'.format(remote_file))

    ## close transport
    transport.close()
