##QV 2019-07-08
##         2019-12-17 integrated in tact
##         2021-02-27 (QV) integrated in acolite renamed from libradtran_run_file

def libradtran_run(runfile):
    import acolite as ac
    import os
    import subprocess
    uvspec = '{}/bin/uvspec'.format(ac.config['libradtran_dir'])

    current_path = os.getcwd()
    binpath = os.path.dirname(uvspec)
    binary = os.path.basename(uvspec)

    os.chdir(binpath)
    outputfile = runfile.replace('.inp', '.out')

    cmd = ['./{}'.format(binary),'< {}'.format(runfile),'> {}'.format(outputfile)]
    cmd = ' '.join(cmd)
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    os.chdir(current_path)
    return(outputfile)
