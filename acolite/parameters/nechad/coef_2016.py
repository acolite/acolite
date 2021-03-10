## def coef_nechad2016
## reads Nechad 2016 calibration data
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-07
## modifications:
##                2018-07-18 (QV) changed acolite import name

def coef_2016():
    import os,sys
    import acolite as ac
    nfile = ac.config['data_dir']+'/Shared/algorithms/Nechad/Nechad_calibration_201609.txt'

    data = []
    with open(nfile, 'r') as f:
        for line in f.readlines():
            if line[0] in [';','!', '/']: continue
            line = line.strip()
            split = line.split('\t')
            if len(split) != 6: continue
            cd = {'par':split[0], 'sensor':split[1], 'band':split[2],
                  'wave':float(split[3]), 'A':float(split[4]), 'C':float(split[5])}
            data.append(cd)

    return(data)
