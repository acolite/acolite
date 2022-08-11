## def bt_lut
## reads ECOSTRESS radiance to BT LUT
## from https://git.earthdata.nasa.gov/projects/LPDUR/repos/ecostress_swath2grid (accessed 2022-08-11)
##
## written by Quinten Vanhellemont, RBINS
## 2022-08-11

def bt_lut(bt_lut_file=None):
    import h5py, os
    import acolite as ac
    lut = None

    ## from https://git.earthdata.nasa.gov/projects/LPDUR/repos/ecostress_swath2grid (accessed 2022-08-11)
    if bt_lut_file is None:
        bt_lut_file = ac.config['data_dir'] + '/ECOSTRESS/EcostressBrightnessTemperatureV01.h5'

    ## read LUT data
    if os.path.exists(bt_lut_file):
        with h5py.File(bt_lut_file, mode='r') as f:
            for b in range(1, 6):
                bk = '{}'.format(b)
                if lut is None: lut = {}
                lut[bk] = f['lut']['radiance_{}'.format(b)][()]
    return(lut)
