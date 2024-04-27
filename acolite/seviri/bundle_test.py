## def bundle_test
## return SEVIRI nat file path from given input bundle
##
## written by Quinten Vanhellemont, RBINS
## 2024-04-12
## modifications:

def bundle_test(bundle_in, match = 'MSG*.nat'):
    import os, glob

    if os.path.isdir(bundle_in):
        bundle_dir = bundle_in
        fn = None
    else:
        bundle_dir = os.path.dirname(bundle_in)
        fn = os.path.basename(bundle_in)

    ## list files
    files = glob.glob('{}/{}'.format(bundle_dir, match))
    files.sort()

    ## ignore xml files for now
    #xfiles = glob.glob('{}/*.xml'.format(bundle_dir))
    #xfiles.sort()

    sel_nat = None
    for file in files:
        if (os.path.splitext(file)[-1] == '.nat'):
            if fn is None:
                fn = os.path.basename(file)
                sel_nat = '{}'.format(file)
            else:
                if fn !=  os.path.basename(file): continue
                sel_nat = '{}'.format(file)

    return(sel_nat)
