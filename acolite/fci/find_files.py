## def find_files
## Find FCI files in given directory
##
## written by Quinten Vanhellemont, RBINS
## 2024-10-18
## modifications: 2025-05-12 (QV) renamed from find_fci_files

def find_files(bundle):
    import glob, os
    import numpy as np
    import acolite as ac

    ## list NetCDFs in inputfile
    files = glob.glob('{}/*.nc'.format(bundle))
    files.sort()

    ## find out which ones are from FCI
    fci_files = {}
    for file in files:
        gatts = ac.shared.nc_gatts(file)
        if (gatts['platform'] in ['MTI1']) & \
           (gatts['data_source'] == 'FCI') & \
           (gatts['processing_level'] == '1C'):
            ## use date_time_position to bundle same observation files
            key = '{}_{}'.format(gatts['date_time_position'], gatts['subtype'])
            if key not in fci_files: fci_files[key] = []
            ## append file and attributes
            fci_files[key].append([file, gatts])

    ## sort on count_in_repeat_cycle
    ## some strips may have processing dates earlier than preceding ones
    for key in fci_files:
        circ = [f[1]['count_in_repeat_cycle'] for f in fci_files[key]]
        idx = np.argsort(circ)
        fci_files[key] = [fci_files[key][i] for i in idx]

    ## return found files per date time position
    return(fci_files)
