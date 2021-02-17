## simple function to import txt config file
## QV 2017-05-24
## Last modifications: 2018-07-25 (QV) changed strip and skip characters

def import_config(file):
    config={}
    with open(file, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if len(line) == 0: continue
            if line[0] in ['#', ';', '%']: continue
            split=line.split('=')
            if len(split) != 2: continue
            split = [s.strip() for s in split]
            config[split[0]]=split[1]
    return(config)
