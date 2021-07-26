## wrapper to launch ACOLITE GUI/CLI
## QV 2018
## last modifications QV 2018-09-12 renamed from acolite.py and added main test
##                    QV 2019-02-21 ignore numpy errors
##                    QV 2021-01-05 added freeze_support call for binary GUI
##                    QV 2021-04-01 updated for generic ACOLITE
##                    QV 2021-05-19 added print of import errors

def launch_acolite():
    ## need to run freeze_support for PyInstaller binary generation
    from multiprocessing import Process, freeze_support
    freeze_support()

    ## import sys to parse arguments
    import sys
    import datetime
    import argparse

    ## import acolite source
    import acolite as ac

    try:
        import acolite as ac
    except:
        print('Could not import ACOLITE source')
        print("Error:", sys.exc_info())
        return()

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
    parser.add_argument('--sensor', help='comma separated sensor list for LUT retrieval', default=None)
    args, unknown = parser.parse_known_args()

    if '--retrieve_luts' in sys.argv:
        ac.acolite.acolite_luts(sensor=args.sensor)
        return()

    ## command line processing, run acolite_run directly
    if '--cli' in sys.argv:
        time_start = datetime.datetime.now() ## time of processing start

        ## add input/output if changed
        inputfile, output = None, None
        if args.inputfile is not None: inputfile = args.inputfile.split(',')
        if args.output is not None: output = args.output

        ## run processing
        if args.settings is None:
            print('No settings file given')
            return()

        ac.acolite.acolite_run(args.settings, inputfile=inputfile, output=output)
    else:
        ret = ac.acolite.acolite_gui(sys.argv, version=ac.version)
        return()

if __name__ == '__main__':
    launch_acolite()
