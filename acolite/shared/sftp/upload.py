## def sftp.upload
## uploads all files in path to sftp, currently only with stat file size check
## written by Quinten Vanhellemont, RBINS
## 2026-04-20
## modifications: 2026-04-20 (QV) added skip paths

def upload(host, local_path, port = 22, remote_path = '/',
            override = False, machine_append = '', skip_paths = [], ):

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

    ## find local files
    local_files = glob.glob('{}/**'.format(local_path), recursive = True)
    local_files.sort()
    print('Recursively found {} files in {}'.format(len(local_files), local_path))
    if len(local_files) == 0: return

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
            else:
                print('File {} exists with correct size'.format(remote_file))

    ## close transport
    transport.close()
