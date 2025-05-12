## def read_tiles
## Reads FCI image tiles
##
## images are assumed to be repacked, or fcidecomp support should be available
##
## written by Quinten Vanhellemont, RBINS
## 2024-10-18
## modifications: 2025-05-12 (QV) moved to function, added dataset support, added flip column and row ranges

def read_tiles(fci_files, dtp, dataset, flip = True, column_range = None, row_range = None):
    import os
    from netCDF4 import Dataset
    import numpy as np

    if 'HRFI' in dtp:
        datatype = 'HRFI'
        datasets = ['vis_06_hr', 'nir_22_hr', 'ir_38_hr', 'ir_105_hr']
        full_dimensions = 22272, 22272
        ssd = 0.5
    elif 'FDHSI' in dtp:
        datatype = 'FDHSI'
        datasets = ['vis_04','vis_05','vis_06','vis_08','vis_09', 'nir_13','nir_16','nir_22',
                    'ir_38','wv_63','wv_73','ir_87','ir_97','ir_105','ir_123','ir_133']
        full_dimensions = 11136, 11136
        ssd = 1.0
    else:
        print('{} not recognised'.format(dtp))
        return

    if dataset not in datasets:
        print('{} not in datasets: {}'.format(dataset, ', '.join(datasets)))
        return

    if column_range is not None:
        if len(column_range) != 2:
            print('Provide two elements for column range')
            return

    if row_range is not None:
        if len(row_range) != 2:
            print('Provide two elements for row range')
            return
            
        ## compute row range in upside down image
        upside_down_row_range = full_dimensions[0]-row_range[1]-1, full_dimensions[0]-row_range[0]-1

    ## run through tiles
    full_row_shape = 0
    data = None
    irr = None
    for t in fci_files[dtp]:
        file, gatts = t
        if '-TRAIL-' in file: continue

        ## read data
        tmp = None
        with Dataset(file) as nc:
            if dataset in nc.groups['data'].groups:
                ## compute shape and range of current row
                curr_row_shape = nc.groups['data'].groups[dataset].groups['measured']['effective_radiance'].shape[0]
                full_row_shape += curr_row_shape
                cur_row_range = full_row_shape-curr_row_shape, full_row_shape

                ## subset current strip
                if row_range is not None:
                    if cur_row_range[1] < upside_down_row_range[0]: continue
                    if cur_row_range[0] > upside_down_row_range[1]: continue

                    ## get y subset in current row
                    y_start = upside_down_row_range[0]-cur_row_range[0]
                    y_end = y_start + upside_down_row_range[1] - upside_down_row_range[0]
                    if y_start < 0: y_start = 0
                    if y_end > curr_row_shape: y_end = curr_row_shape

                    print('Reading {} strip {} {}:{}'.format(dataset, gatts['count_in_repeat_cycle'], y_start, y_end))
                    if column_range is not None:
                        tmp = nc.groups['data'].groups[dataset].groups['measured']['effective_radiance'][y_start:y_end, column_range[0]:column_range[1]]
                    else:
                        tmp = nc.groups['data'].groups[dataset].groups['measured']['effective_radiance'][y_start:y_end, :]
                else:
                    if column_range is not None:
                        tmp = nc.groups['data'].groups[dataset].groups['measured']['effective_radiance'][:, column_range[0]:column_range[1]]
                        print(tmp.shape)
                    else:
                        tmp = nc.groups['data'].groups[dataset].groups['measured']['effective_radiance'][:]

                ## read irradiance
                if irr is None:
                    irr = nc.groups['data'].groups[dataset].groups['measured']['channel_effective_solar_irradiance'][:]
                else:
                    if irr != nc.groups['data'].groups[dataset].groups['measured']['channel_effective_solar_irradiance'][:]:
                        print('Irradiance does not match')

                tmp[tmp.mask] = np.nan
                tmp = tmp.astype(np.float32)

        ## continue if missing dataset
        if tmp is None:
            print('Dataset {} not in {}'.format(dataset, file))
            continue

        ## create disk
        if data is None:
            data = tmp
            line = np.zeros(tmp.shape)
        else:
            data = np.vstack((data, tmp))

    ## flip disk North side up
    if (flip) & (data is not None): data = np.flipud(data)

    return(data, irr)
