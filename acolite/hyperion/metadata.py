## def metadata
## read Hyperion metadata
##
## written by Quinten Vanhellemont, RBINS
## 2018-01-30
## modifications: 2021-08-04 (QV) added metadata file search, renamed and added to generic acolite

def metadata(bundle):
    import glob
    mdata = {}
    tmp = glob.glob('{}/*MTL*.TXT'.format(bundle))
    if len(tmp) == 1:
        metafile = tmp[0]
        with open(metafile,'r') as f:
            for l,line in enumerate(f.readlines()):
                if line[0:3]=='END': continue
                tag, value = line.split('=')
                tag = tag.strip()
                value = value.strip()
                if '"' in value: value = value.strip('"')

                if tag == 'GROUP':
                    group = {'name':value}
                elif tag == 'END_GROUP':
                    mdata[group['name']]=group
                else:
                    group[tag]=value
    return(mdata)
