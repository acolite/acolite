## def acolite_luts
## function to retrieve & test required LUTs, e.g. to set up Docker
## written by Quinten Vanhellemont, RBINS
## 2021-07-26
## modifications:
##

def acolite_luts(sensor = None, hyper = False, pars = ['romix', 'romix+rsky_t']):
    import acolite as ac

    ## get all RSR if no sensor(s) specified
    if (sensor is None) | (sensor.upper() == 'NONE'):
        rsrd = ac.shared.rsr_dict()
    else:
        if type(sensor) is not list:
            if ',' in sensor:
                sensor = [s.strip() for s in sensor.split(',')]
            elif ' ' in sensor:
                sensor = [s for s in sensor.split(' ') if len(s) > 0]
            else:
                sensor = [sensor]
        rsrd = {}
        for s in sensor:
            s = s.upper()
            if s in ['CHRIS', 'PRISMA', 'HYPER']:
                hyper = True
                continue

            rd = ac.shared.rsr_dict(sensor=s)
            if len(rd) == 0:
                print('Sensor {} not recognised.'.format(s))
                continue
            rsrd[s] = rd[s]

    sensors = list(rsrd.keys())
    sensors.sort()

    ## Add "None" to sensors to retrieve generic LUT
    if hyper: sensors += [None]

    ## run through wanted sensors
    for s in sensors:
        if s in ['L5_TM_B6', 'L7_ETM_B6', 'L8_OLI_ORANGE', 'L8_TIRS']: continue ## skip thermals

        print('Testing {}'.format('sensor {}'.format(s) if s is not None else 'generic LUT'))

        ## try getting gas transmittance
        tg_dict = ac.ac.gas_transmittance(0, 0, uoz=0.3, uwv=1.6, rsr=None if s is None else rsrd[s]['rsr'])

        ## get sensor LUT
        tmp = ac.aerlut.import_luts(sensor = s)

        ## get reverse LUT
        if (s is not None) & (s in ['L5_TM', 'L7_ETM', 'L8_OLI', 'S2A_MSI', 'S2B_MSI', 'S3A_OLCI', 'S3B_OLCI']):
            for par in pars:
                revl = ac.aerlut.reverse_lut(s, par=par)
