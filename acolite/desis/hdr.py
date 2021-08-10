## Read DESIS ENVI hdr file
## written by Quinten Vanhellemont, RBINS
## 2021-08-10

def hdr(headerfile):
    with open(headerfile, 'r', encoding='utf-8') as f:
        tag = None
        header = {}
        for il, line in enumerate(f.readlines()):
            line = line.strip()
            sp = [s.strip().strip('{},') for s in line.split('=')]

            if len(sp) == 1:
                if tag is None: continue
                if len(sp[0]) == 0: continue
                header[tag].append(sp[0])
            else:
                tag = sp[0]
                if tag not in header: header[tag] = []
                if len(sp[1]) == 0: continue
                if len(sp)>2:
                    header[tag].append(sp[1:])
                else:
                    header[tag].append(sp[1])

    for k in header:
        if len(header[k]) == 1: header[k] = header[k][0]
        if k in ['samples', 'lines', 'bands', 'header offset',
                 'data type', 'byte order', 'data ignore value']:
            header[k] = int(header[k])
        if k in ['wavelength', 'fwhm', 'data gain values', 'data offset values']:
            header[k] = [float(v) for v in header[k]]
    return(header)
