## def coef_chl_oc
## reads chl_oc calibration data
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-07
## modifications:
##                2018-07-18 (QV) changed acolite import name

def coef(config='defaults'):
    import os,sys
    import acolite as ac
    nfile = ac.config['data_dir']+'/Shared/algorithms/chl_oc/{}.txt'.format(config)

    data = {}
    with open(nfile, 'r') as f:
        for line in f.readlines():
            if line[0] in [';','!', '/']: continue
            line = line.strip()
            split = line.split('\t')
            if len(split) != 6: continue
            cd = {'par':split[0],
                  'blue':[float(i) for i in split[1].split(',')],
                  'green':[float(i) for i in split[2].split(',')],
                  'chl_coef':[float(i) for i in split[3].split(',')], 'reference':(split[4]), 'sensor':split[5]}
            if cd['sensor'] not in data: data[cd['sensor']] = {}
            data[cd['sensor']][cd['par']] = cd
    return(data)
