## def acolite_luts
## function to retrieve & test required LUTs, e.g. to set up Docker
## written by Quinten Vanhellemont, RBINS
## 2021-07-26
## modifications:
##                2021-10-24 (QV) added LUT identifiers and pressures, get_remote keyword
##                2022-04-12 (QV) add par parameter to import_luts

def acolite_luts(sensor = None, hyper = False,
                 get_remote = True, compute_reverse = True,
                 pressures = [500, 750, 1013, 1100],
                 base_luts = ['ACOLITE-LUT-202110-MOD1', 'ACOLITE-LUT-202110-MOD2'],
                 rsky_lut = 'ACOLITE-RSKY-202102-82W',
                 pars = ['romix', 'romix+rsky_t']):
    import acolite as ac

    ## get all RSR if no sensor(s) specified
    if type(sensor) is str:
        if sensor.upper() == 'NONE': sensor = None

    if (sensor is None):
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
            if (s.upper() in ['CHRIS', 'PRISMA', 'HYPER']) or ('DESIS' in s.upper()):
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
        if s is not None:
            if s in ['L5_TM_B6', 'L7_ETM_B6', 'EO1_ALI_ORANGE',
                     'L8_OLI_ORANGE', 'L8_TIRS', 'L9_OLI_ORANGE', 'L9_TIRS']: continue ## skip thermals
            if 'DESIS' in s: continue ## skip DESIS

        print('Testing {}'.format('sensor {}'.format(s) if s is not None else 'generic LUT'))

        ## try getting gas transmittance
        tg_dict = ac.ac.gas_transmittance(0, 0, uoz=0.3, uwv=1.6, rsr=None if s is None else rsrd[s]['rsr'])

        ## get sensor LUT
        tmp = ac.aerlut.import_luts(sensor = s, get_remote = get_remote, pressures = pressures, par=pars[-1],
                                    base_luts = base_luts, rsky_lut = rsky_lut)

        ## get reverse LUT
        if (compute_reverse) & (s is not None) & (s in ['L5_TM', 'L7_ETM', 'L8_OLI', 'L9_OLI', 'S2A_MSI', 'S2B_MSI',
                                                        'S3A_OLCI', 'S3B_OLCI', 'EN1_MERIS']):
            for par in pars:
                revl = ac.aerlut.reverse_lut(s, get_remote = get_remote, par=par, pressures = pressures,
                                            base_luts = base_luts, rsky_lut = rsky_lut)
