## def metadata_read
## reads landsat metadata
## written by Quinten Vanhellemont, RBINS
## 2017-04-13
## modifications: 2018-09-19 QV changed strip, added encoding

def metadata_read(metafile):
    import acolite as ac

    mdata={}
    with open(metafile,'r', encoding="utf-8") as f:
        for line in f.readlines():
            line = line.strip()
            split = line.split('=')
            for s, sp in enumerate(split):
                split[s] = sp.strip('" ')
            if len(split) != 2: continue
            if split[0] == 'GROUP':
                cur_group = split[1]
                group_data = {}
            elif split[0] == 'END_GROUP':
                mdata[cur_group] = group_data
            else:
                group_data[split[0]] = split[1]

    ## do updates for ALI
    mdata = ac.landsat.metadata_ali(mdata)

    return(mdata)
