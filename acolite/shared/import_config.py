## simple function to import txt config file
## QV 2017-05-24
## Last modifications: 2018-07-25 (QV) changed strip and skip characters
##                     2023-07-30 (QV) added parse keyword

def import_config(file, parse=False):
    import os
    config={}
    if os.path.exists(file):
        with open(file, 'r', encoding='utf-8') as f:
            for line in f.readlines():
                line = line.strip()
                if len(line) == 0: continue
                if line[0] in ['#', ';', '%']: continue
                first = line.find('=')
                if first <1: continue
                split = line[0:first], line[first+1:]
                split = [s.strip() for s in split]
                if parse:
                    split[1] = split[1].split(',')
                    try:
                        config[split[0]] = [float(s) for s in split[1]]
                    except:
                        config[split[0]] = split[1]
                    if len(config[split[0]]) == 1: config[split[0]]=config[split[0]][0]
                    if config[split[0]] in ['True','true']: config[split[0]]=True
                    if config[split[0]] in ['False','false']: config[split[0]]=False
                    if config[split[0]] in ['None','none']: config[split[0]]=None
                else:
                    config[split[0]]=split[1]
    else:
        print('{} not found'.format(file))
    return(config)
