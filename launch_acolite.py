## wrapper to launch ACOLITE GUI/CLI
## QV 2018
## last modifications QV 2018-09-12 renamed from acolite.py and added main test
##                    QV 2019-02-21 ignore numpy errors
##                    QV 2021-01-05 added freeze_support call for binary GUI
##                    QV 2021-04-01 updated for generic ACOLITE

def launch_acolite():
    ## need to run freeze_support for PyInstaller binary generation
    from multiprocessing import Process, freeze_support
    freeze_support()

    ## import acolite source
    try:
        import acolite as ac
    except:
        print('Could not import ACOLITE source')
        return()

    ## import sys to parse arguments
    import sys
    import datetime
    import argparse

    ## fix matplotlib backend to Agg
    ## skip import if --nogfx is given
    if not (('--cli' in sys.argv) & ('--nogfx' in sys.argv)):
        import matplotlib
        matplotlib.use("Agg")

    ## ignore numpy errors
    import numpy as np
    olderr = np.seterr(all='ignore')

    ## run command line if --cli provided, otherwise use gui
    parser = argparse.ArgumentParser(description='ACOLITE')
    parser.add_argument('--settings', help='settings file', default=None)
    parser.add_argument('--inputfile', help='list of images', default=None)
    parser.add_argument('--output', help='output directory', default=None)
    args, unknown = parser.parse_known_args()

    ## command line processing, run acolite_run directly
    if '--cli' in sys.argv:
        time_start = datetime.datetime.now() ## time of processing start
        
        ## add input/output if changed
        inputfile, output = None, None
        if args.inputfile is not None: inputfile = args.inputfile.split(',')
        if args.output is not None: output = args.output

        ## run processing
        ac.acolite.acolite_run(args.settings, inputfile=inputfile, output=output)
    else:
        ret = ac.acolite.acolite_gui(sys.argv, version=ac.version)
        return()

if __name__ == '__main__':
    launch_acolite()
