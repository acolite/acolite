## hydrolight.config.write
## write HE5 config to Iroot.txt file
##
## written by Quinten Vanhellemont, RBINS
## 2018-04-23
## modifications: 2026-06-16 (QV) new function

def write(file, config):
    import os
    fdir = os.path.dirname(file)
    if not os.path.exists(fdir): os.makedirs(fdir)

    with open(file, 'w') as f:
        for i in range(0,13):
            rec = i+1
            reclen = 1
            lines = []

            if rec in [1,2,3,4,7,9,10,12]:
                if rec == 1: recv = [["icompile", "Parmin", "Parmax", "PhiChl",
                                      "Raman0", "RamanXS", "iDynZ", "RamanExp" ]]
                if rec == 2: recv = [["ititle"]]
                if rec == 3: recv = [["rootname"]]
                if rec == 4:
                    recv = [["iOptPrnt","iOptDigital","iOptExcelS",
                             "iOptExcelM","iOptRad","nwskip"],
                            ["iIOPmodel","iSkyRadModel","iSkyIrradModel","iChl","iCDOM"]]
                if rec == 7:
                    recv = [["ibiolum","ichlfl","icdomfl","iraman","icompchl"]]
                if rec == 9:
                    recv = [["windspd","refr","temp","salinty"]]
                if rec == 10:
                    recv = [["ibotm","rflbot"]]
                if rec == 12:
                    recv = [['PureWaterDataFile'],['nac9Files'],
                            ['ac9DataFile'],['Ac9FilteredDataFile'],
                            ['HydroScatDataFile'],['ChlzDataFile'],
                            ['CDOMDataFile'],['RbottomFile'],['TxtDataFile'],
                            ['IrradDataFile'],['S0biolumFile']]
                for rv in recv:
                    lines.append(['{}'.format(config[v]['value']) for v in rv])
            elif rec == 5:
                rv = ['ncomp','nconc']
                lines.append(['{}'.format(config[v]['value']) for v in rv])
                lines.append(config['compconc'])

                for i in range(0,config['nconc']['value']):
                    lines.append(config['abspar{}'.format(i+1)])
                for i in range(0,config['nconc']['value']):
                    lines.append(config['absfile{}'.format(i+1)])

                for i in range(0,config['ncomp']['value']):
                    lines.append(config['scapar{}'.format(i+1)])
                for i in range(0,config['ncomp']['value']):
                    lines.append(config['scafile{}'.format(i+1)])

                for i in range(0,config['ncomp']['value']):
                    lines.append(config['phapar{}'.format(i+1)])
                for i in range(0,config['ncomp']['value']):
                    lines.append(config['phafile{}'.format(i+1)])
            elif rec == 6:
                from math import ceil
                rv = ['nwave']
                lines.append(['{}'.format(config[v]['value']) for v in rv])
                nw = config['nwave']['value']
                for i in range(0, ceil(nw/12)):
                    lines.append(config['waverec{}'.format(i+1)])
            elif rec == 8:
                lines.append(config['skyrec{}'.format(1)])
                rv=["jday","rlat","rlon", "pres",
                      "am","rh","wv","vi","wsm","ro3"]
                lines.append(['{}'.format(config[v]['value']) for v in rv])
            elif rec == 11:
                dtags = [dt for dt in config.keys() if 'depthrec' in dt]
                dtags.sort()
                for dt in dtags: lines.append(config[dt])
                #lines.append(config['depthrec{}'.format(1)])
            else: continue

            for line in lines:
                if type(line) == list: line = ', '.join(line)
                f.write(line)
                f.write('\n')
