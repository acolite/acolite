# def libradtran_run
## simple script to run libRadtran uvspec
##
## written by Quinten Vanhellemont, RBINS
## 2019-07-08
## modifications: 2019-12-17 (QV) integrated in tact
##                2021-02-27 (QV) integrated in acolite renamed from libradtran_run_file
##                2025-05-20 (QV) check if system uvspec binary is available

def libradtran_run(runfile):
    import acolite as ac
    import os, shutil
    import subprocess

    outputfile = runfile.replace('.inp', '.out')

    current_path = None
    if ac.settings['run']['use_system_libradtran']:
        uvspec = shutil.which('uvspec')
    else:
        uvspec = None

    ## use uvspec in libradtran_dir/bin
    if uvspec is None:
        uvspec = '{}/bin/uvspec'.format(ac.config['libradtran_dir'])
        current_path = os.getcwd()
        binpath = os.path.dirname(uvspec)
        binary = os.path.basename(uvspec)
        os.chdir(binpath)
        cmd = ['./{}'.format(binary),'< "{}"'.format(runfile),'> "{}"'.format(outputfile)]
    else:
        cmd = ['{}'.format(uvspec),'< "{}"'.format(runfile),'> "{}"'.format(outputfile)]

    ## join command
    cmd = ' '.join(cmd)
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)

    ## change back to current path
    if os.getcwd() != current_path: os.chdir(current_path)
    return(outputfile)
