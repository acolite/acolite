## QV Jul 2019
##         2019-12-17 renamed, integrated in tact
##         2021-02-27 (QV) integrated in acolite renamed from cfg_create
##         2022-08-11 (QV) extended wavelength range, added option to use reptran if available

def libradtran_cfg(thv=0.0, albedo=0.0, look_down=True,
               ths=0.0, phi=0.0, phi0=0.0, atmosphere=None, wavelength = [7000, 14000],
               reptran = 'medium', zout = ['SUR', 'TOA'], parameters = ['lambda','eup','uu'],
               radiosonde = None, sur_temperature=None, brightness=False, runfile=None):

        import os
        import numpy as np
        import acolite as ac

        if type(zout) is str: zout = [zout]

        config=[
            "data_files_path {}/data".format(ac.config['libradtran_dir']),
            'source thermal',
            'wavelength {} {}'.format(wavelength[0], wavelength[1]),
            'albedo {}'.format(albedo),
            'rte_solver disort',
            'zout {}'.format(' '.join(zout)),
            'sza {}'.format(ths),
            'phi {}'.format(phi),
            'phi0 {}'.format(phi0),
            'umu {}'.format(np.cos(thv*(np.pi/180.))*(1 if look_down else -1)),
            'output_user {}'.format(' '.join(parameters)),
            'quiet'
        ]

        ## add reptran if exists
        repfine = ac.config['libradtran_dir'] + '/data/correlated_k/reptran/reptran_solar_fine.cdf'
        if os.path.exists(repfine): config.insert(2, 'mol_abs_param reptran {}'.format(reptran))

        if sur_temperature is not None: config += ['sur_temperature {}'.format(sur_temperature)]
        if atmosphere is not None: config +=['atmosphere_file {}'.format(atmosphere)]
        if brightness: config += ['output_quantity brightness']
        else: config += ['output_process per_nm']
        if radiosonde is not None:
            config += ['radiosonde {} H2O RH'.format(radiosonde)]

        if runfile is not None:
            with open(runfile, 'w') as file:
                for line in config: file.write(line+'\n')
        return(config)
