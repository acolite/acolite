## def settings_read
## read ACOLITE settings from file to "settings" dict
##
## Written by Quinten Vanhellemont 2017-11-30
## Last modifications: 2018-06-06 QV added None check for limit
##                     2018-09-12 QV added length check for splitting on =
##                     2018-09-19 QV added encoding
##                     2020-07-14 QV added strip to parameters

def read(file):
    import os
    settings={}
    if os.path.exists(file):
        with open(file, 'r', encoding="utf-8") as f:
            for line in f.readlines():
                line = line.strip()
                if len(line) == 0: continue
                if line[0] in ['#',';']: continue
                split = line.split('=')
                if len(split) < 2: continue
                split = [s.strip() for s in split]
                var = split[0]
                val = [s.strip() for s in split[1].split(',')]
                if len(val) == 1:
                    val = val[0]
                    if val in ['True','true','1']: val=True
                    if val in ['False','false','0']: val=False
                    if val in ['None','none']: val=None
                if (var in ['limit']) & (val is not None): val = [float(i) for i in val]
                settings[var]=val
    return(settings)
