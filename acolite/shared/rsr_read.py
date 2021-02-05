## def rsr_read
## imports rsr file
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07-11
## modifications: 2016-11-23 QV added support for the S2A RSR file
##                2018-03-27 QV Fixed bug in L8 RSR file on windows
##                2018-09-19 (QV) added encoding

def rsr_read(file=None):
    import fnmatch
    if file is not None:
        with open(file, 'r', encoding='utf-8') as f:
            rwave=[]
            rresp=[]
            bands=[]
            bid = 0
            data = dict()
            for i, line in enumerate(f.readlines()):
                if (line[0][:1] == '#') | (line[0][:1] == ';'):
                    if '/' in line: continue
                    if (fnmatch.fnmatch(line, '*Band*')) | (fnmatch.fnmatch(line, '*BAND*')): 
                        tmp = line.split()
                        band = tmp[-1]
                        if bid > 0:
                            bdata = {'wave':rwave, 'response':rresp}
                            data[prev_band]=bdata
                            bands.append(prev_band)
                        prev_band = band
                        bid+=1
                        rwave=[]
                        rresp=[]

                else:
                    ls = line.split()
                    if float(ls[0]) > 100:
                        rwave.append(float(ls[0])/1000.)
                    else:
                        rwave.append(float(ls[0]))
                    rresp.append(float(ls[1]))

            if len(rwave) != 0:
                bdata = {'wave':rwave, 'response':rresp}
                data[prev_band]=bdata
                bands.append(prev_band)

        return data,bands
