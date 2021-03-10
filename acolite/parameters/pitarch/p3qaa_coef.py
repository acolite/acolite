## function to read Jaime's 3 band QAA coeffs
## QV 2021-02-15

def p3qaa_coef():
    import acolite as ac

    pars = ['alpha','aw','bbw','beta1','beta2','bg_ratio','center_wl','chi','coef_z_sd']

    p3qaa = {}

    for par in pars:
        cfg_file = ac.config['data_dir']+'/Shared/algorithms/P3QAA/{}.csv'.format(par)
        with open(cfg_file, 'r') as f:
            header = None
            for il, line in enumerate(f.readlines()):
                line = line.strip()
                if len(line) == 0: continue
                if line[0] in ['#', ';']: continue
                sp = line.split(',')
                if header is None:
                    header = sp
                else:
                    sensor = sp[0]
                    if sensor not in p3qaa:
                        p3qaa[sensor] = {}
                    if par not in p3qaa[sensor]:
                        p3qaa[sensor][par] = {h:sp[ih] for ih, h in enumerate(header)}

                        for k in ['p4', 'p3','p2','p1', 'p0', 'R', 'G', 'B']:
                            if k in p3qaa[sensor][par]:
                                p3qaa[sensor][par][k] = float(p3qaa[sensor][par][k])
    return(p3qaa)
