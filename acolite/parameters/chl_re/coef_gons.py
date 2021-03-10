## def coef_chl_re_gons
## reads chl_re coefficients for Gons algo
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-27
## modifications:
##                2018-07-18 (QV) changed acolite import name

def coef_gons(config='gons'):
    import os,sys
    import acolite as ac
    nfile = ac.config['data_dir']+'/Shared/algorithms/chl_re/{}.txt'.format(config)

    header = ['algorithm','red_band','rededge_band','nir_band',
              'astar_chl','chl_coef','validity','reference']

    data = {}
    with open(nfile, 'r') as f:
        for line in f.readlines():
            if line[0] in [';','!', '/']: continue
            line = line.strip()
            split = line.split('\t')
            if len(split) != 8: continue
            cd = {}
            for ih, h in enumerate(header):
                if h in ['red_band','rededge_band','nir_band','astar_chl']:
                    cv = float(split[ih])
                elif h in ['chl_coef','validity']:
                    cv = [float(j) for j in split[ih].split(',')]
                else:
                    cv = split[ih]
                cd[h] = cv
            data[cd['algorithm']] = cd
    return(data)
