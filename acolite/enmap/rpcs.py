## def rpcs
## parse EnMAP RPC data
## written by Quinten Vanhellemont, RBINS
## 2022-09-21
## modifications: 2022-09-22 (QV) added function

def rpcs(band_data, bname='1'):
    keys = {'HEIGHT_OFF': 'HEIGHT_OFF',
        'HEIGHT_SCALE': 'HEIGHT_SCALE',
        'LAT_OFF': 'LAT_OFF',
        'LAT_SCALE':'LAT_SCALE',
        'LINE_DEN_COEFF': 'ROW_DEN_*',
        'LINE_NUM_COEFF': 'ROW_NUM_*',
        'LINE_OFF': 'ROW_OFF',
        'LINE_SCALE': 'ROW_SCALE',
        'LONG_OFF': 'LONG_OFF',
        'LONG_SCALE': 'LONG_SCALE',
        'SAMP_DEN_COEFF': 'COL_DEN_*',
        'SAMP_NUM_COEFF': 'COL_NUM_*',
        'SAMP_OFF': 'COL_OFF',
        'SAMP_SCALE': 'COL_SCALE',
       }

    RPCs = {}
    for k in keys:
        if keys[k][-1] == '*':
            cur = keys[k][0:-1]
            pars = [str(band_data[bname][c]) for c in band_data[bname] if c[0:len(cur)] == cur]
            RPCs[k] = ' '.join(pars)
        else:
            RPCs[k] = str(band_data[bname][keys[k]])

    return(RPCs)
