## he5.file.run
## runs Hydrolight for given config file
##
## written by Quinten Vanhellemont, RBINS
## 2018-04-23
## modifications: 2022-11-26 (QV) added ECOLIGHT and remove keyword
##                2026-06-16 (QV) new function, removed copy of input as now simulations are copied to the input unique run directory
##                2026-06-18 (QV) use hl dir from config

def run(runfile, result_dir = None,
            he5_dir = None, he5_bin = "hlbin", force = False, monitor = False,
            copy_printout = False, copy_excel = False, remove = False):
    import uuid, os, subprocess
    from shutil import copy
    import acolite as ac

    ## get he5_dir from config
    if he5_dir is None: he5_dir = ac.config['hydrolight_dir']

    ## get current directory and make unique subdirectory
    #cwd = os.getcwd()
    #sbd = '.{}'.format(uuid.uuid4())
    #os.mkdir(sbd)
    cwd = os.getcwd()
    os.chdir("{}/{}".format(he5_dir, 'Code'))

    lockfile = ".he5_running"
    if he5_bin=="hlbin":
        hn = 'Hydrolight'
    if he5_bin=="elbin":
        hn = 'Ecolight'

    if (not os.path.exists(lockfile)) | (force):
        # make lockfile
        if not force:
            with open(lockfile, 'w') as f:
                f.write("")

        ## run sixs command
        if not monitor:
            cmd = './{} < {}'.format(he5_bin, runfile)
            sp = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)

        else:
            #cmd = ['./{} < '.format(he5_bin), runfile]
            #sp = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            cmd = './{} < {}'.format(he5_bin, runfile)
            sp = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)
            print(sp.stdout.decode())
            # Grab stdout line by line as it becomes available.  This will loop until
            # p terminates.
            #while sp.poll() is None:
            #    l = sp.stdout.readline() # This blocks until it receives a newline.
            #    print(l.strip())
            # When the subprocess terminates there might be unconsumed output
            # that still needs to be processed.
            #print(sp.stdout.read())

        # delete lockfile
        if os.path.exists(lockfile): os.remove(lockfile)

        # copy output file
        root = os.path.splitext(os.path.basename(runfile))[0][1:]
        dresult = '{}/{}/{}/{}/{}'.format(he5_dir, 'output', hn, 'digital', 'D{}.txt'.format(root))
        if os.path.exists(dresult):
            if result_dir is not None:
                if not os.path.exists(result_dir): os.makedirs(result_dir)
                copy(dresult, result_dir)
                #copy(runfile, result_dir)
                if remove: os.remove(dresult)

                if copy_printout:
                    presult = '{}/{}/{}/{}/{}'.format(he5_dir, 'output', hn, 'printout', 'P{}.txt'.format(root))
                    if os.path.exists(presult):
                        copy(presult, result_dir)
                        if remove: os.remove(presult)

                if copy_excel:
                    eresult = '{}/{}/{}/{}/{}'.format(he5_dir, 'output', hn, 'excel', 'S{}.txt'.format(root))
                    if os.path.exists(eresult):
                        copy(eresult, result_dir)
                        if remove: os.remove(eresult)

                    eresult = '{}/{}/{}/{}/{}'.format(he5_dir, 'output', hn, 'excel', 'M{}.txt'.format(root))
                    if os.path.exists(eresult):
                        copy(eresult, result_dir)
                        if remove: os.remove(eresult)

                if True:
                    lresult = '{}/{}/{}/{}/{}'.format(he5_dir, 'output', hn, 'digital', 'L{}.txt'.format(root))
                    if os.path.exists(lresult):
                        #copy(lresult, result_dir)
                        if remove: os.remove(lresult)

    ## go back to current directory and remove unique simulation directory
    os.chdir(cwd)
