## def optimise_read
## read target data for optimise functions
## written by Quinten Vanhellemont, RBINS
## 2025-03-17
## modifications: 2025-03-17 (QV) split off from optimise_aot_homogeneous

def optimise_read(setu):
    import acolite as ac
    import numpy as np

    target_data_read = False
    ## if optimise_target_rhos_file is given, read the data, and resample to the sensor RSR
    if (setu['optimise_target_rhos_file'] is not None):
        print('Reading target rhos from file: {}'.format(setu['optimise_target_rhos_file']))
        if (setu['optimise_target_rhos_file_type'].lower() in ['csv', 'tsv']):
            delimiter = {'csv': ',', 'tsv': '\t'}[setu['optimise_target_rhos_file_type']]
            data_import = np.loadtxt(setu['optimise_target_rhos_file'], delimiter = delimiter, dtype = np.float32)
            if len(data_import.shape) != 2:
                print('Wrong shape of data in {}'.format(setu['optimise_target_rhos_file']))
                print(data_import.shape)
            ## check dimensions of data
            if data_import.shape[0] == 2:
                wave_data = data_import[0,:]
                rhos_data = data_import[1,:]
                target_data_read = True
            elif data_import.shape[1] == 2:
                wave_data = data_import[:,0]
                rhos_data = data_import[:,1]
                target_data_read = True
            else:
                print('Wrong shape of data in {}'.format())
                print(data_import.shape)
        elif (setu['optimise_target_rhos_file_type'].lower() in ['nc', 'netcdf']):
            tg, td, ta = ac.shared.nc_read_all(setu['optimise_target_rhos_file'])
            if (setu['optimise_target_rhos_file_wavelength'] not in td):
                print('Data type optimise_target_rhos_file_wavelength={} not in {}'.format(setu['optimise_target_rhos_file_wavelength'], setu['optimise_target_rhos_file']))
            elif (setu['optimise_target_rhos_file_reflectance'] not in td):
                print('Data type optimise_target_rhos_file_reflectance={} not in {}'.format(setu['optimise_target_rhos_file_reflectance'], setu['optimise_target_rhos_file']))
            else:
                wave_data = td[setu['optimise_target_rhos_file_wavelength']].flatten()
                rhos_data = td[setu['optimise_target_rhos_file_reflectance']].flatten()
                target_data_read = True
                del tg, td, ta
        else:
            print('Data type optimise_target_rhos_file_type={} not configured'.format(setu['optimise_target_rhos_file_type']))

    ## return data
    if not target_data_read:
        return
    else:
        return(wave_data, rhos_data)
