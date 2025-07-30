## def get_radiometric_calibration
## read radiometric calibration for Huanjing
##
## written by Quinten Vanhellemont, RBINS
## 2025-07-30
## modifications:

def get_radiometric_calibration(calibration_file=None):
    import acolite as ac
    if calibration_file is None:
        calibration_file = ac.config['data_dir'] + '/Huanjing/HJ2_radiometric_calibration.txt'
    radiometric_calibration = {}
    with open(calibration_file, 'r', encoding = 'utf-8') as f:
        for il, line in enumerate(f.readlines()):
            line = line.strip()
            sp = line.split(',')
            if sp[0] not in radiometric_calibration: radiometric_calibration[sp[0]] = {}
            radiometric_calibration[sp[0]][sp[1]] = [float(v) for v in sp[2:]]
    return(radiometric_calibration)
