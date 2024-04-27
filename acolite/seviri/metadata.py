## def metadata
## get metadata from SEVIRI nat file
##
## written by Quinten Vanhellemont, RBINS
## 2024-04-09
## modifications:

def metadata(file, quiet = True):
    meta = {}
    with open(file, 'rb') as f:
        read = []
        for il, byteline in enumerate(f.readlines()):
            if il == 48: break
            read.append(byteline)
            line = byteline.decode().replace('\x00', '').strip() ## decode bytestring
            if not quiet: print(line)
            sp = [s.strip() for s in line.split(':')]
            if ' ' in sp[1]: sp[1] = ' '.join(sp[1].split()) ## remove multiple spaces
            meta[sp[0]] = sp[1]
    return(meta)
