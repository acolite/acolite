## read CHRIS interband calibration data
## source file INTERBAND_CALIBRATION_COEFFS_CHRIS_PROBA.csv from Héloïse Lavigne
##
## QV 2021-05-21

def interband_calibration(file=None):
    import re
    import dateutil.parser, datetime
    import acolite as ac

    if file is None:
        file = '{}/CHRIS/INTERBAND_CALIBRATION_COEFFS_CHRIS_PROBA.csv'.format(ac.config['data_dir'])

    ## read data
    header = None
    periods = {}
    with open(file, 'r') as f:
        for il, line in enumerate(f.readlines()):
            line = line.strip()

            if line[0:9] == '# Period ':
                period = line.split()[2]
                dates = line.split(':')[-1].strip().split(' to ')

                ## get year from second date if not given
                year = None
                for i in range(2):
                    for di, d in enumerate(dates):
                        match = re.match(r'.*([1-3][0-9]{3})', d)
                        if match is not None:
                            year = match.group(1)
                        else:
                            if year is None: continue
                            dates[di] = '{} {}'.format(d, year)

                ## parse dates
                dr = [dateutil.parser.parse(d) for d in dates]

                ## set day to first and last of month
                dt0 = datetime.datetime(dr[0].year, dr[0].month, 1)
                ml1 = (dt0.replace(month = dt0.month % 12 + 1, day = 1)-datetime.timedelta(days=1)).day
                dt1 = datetime.datetime(dr[1].year, dr[1].month, 1)+datetime.timedelta(days=ml1)
                periods[period] = [dt0, dt1]
            elif line[0] in ['#']:
                continue
            else:
                sp = line.split(',')
                if header is None:
                    header = [s for s in sp]
                    data = {}
                else:
                    v = {h:sp[ih] for ih,h in enumerate(header)}
                    for h in v:
                        if h not in ['band']: v[h] = float(v[h])
                    data[v['band']] = v

    ## convert to dict
    datap = {}
    for b in data:
        for p in periods:
            if p not in datap: datap[p] = {'data':{}, 'range': periods[p]}
            datap[p]['data'][b] = {'cal':data[b]['CAL_Period_{}'.format(p)], 'std':'SD_Period_{}'.format(p)}
    return(datap)
