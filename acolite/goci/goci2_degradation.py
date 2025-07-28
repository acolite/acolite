## def goci2_degradation
## reads GOCI2 degradatation info
##
## written by Quinten Vanhellemont, RBINS
## 2025-07-24
## modifications:

def goci2_degradation(date, t0 = '2020-03-10T00:00:00Z'):
    import numpy as np
    import acolite as ac
    import datetime, dateutil.parser

    file = ac.config['data_dir']+'/GOCI/GOCI2_sensor_degradation.csv'
    data = np.loadtxt(file, delimiter = ',')

    datetime_t0 = dateutil.parser.parse(t0)
    timedelta = dateutil.parser.parse(date) - datetime_t0
    days_since_t0 = timedelta.days + timedelta.seconds / 3600

    if days_since_t0 < data[0,0]:
        print('Warning: {} before T0 in {}'.format(date, file))
    if days_since_t0 > data[-1,0]:
        print('Warning: {} after last time since T0 in {}'.format(date, file))

    coefficients = []
    for b in range(12):
        cur = np.interp(days_since_t0, data[:,0], data[:, b+1])
        coefficients.append(cur)

    return(np.asarray(coefficients))
