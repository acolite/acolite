## hydrolight.config.conf4c
## HE5.2 input configuration dictionary based on HE52TechDoc
## with QV defaults
##
## written by Quinten Vanhellemont, RBINS
## 2018-04-23
## modifications: 2021-09-29 (QV) new 4 component version changed defaults/defaults handling
##                2026-06-16 (QV) new function, changed wave defaults from step = 10, min_wave = 400, max_wave = 1500

def conf4c(step = 5, min_wave = 400, max_wave = 900, **kwargs):
        import numpy as np

        rdefs = {## record 1 -
                 "icompile":{'value':0, 'option':''},
                 "Parmin":{'value':400, 'option':''},
                 "Parmax":{'value':700, 'option':''},
                 "PhiChl":{'value':0.02, 'option':''},
                 "Raman0":{'value':488, 'option':''},
                 "RamanXS":{'value':0.00026, 'option':''},
                 "iDynZ":{'value':1, 'option':''},
                 "RamanExp":{'value':5.3, 'option':''},
                 ## record 2 title
                 "ititle":{'value':'HydroLight 5.3 Input File', 'option':''},
                 ## record 3 rootname
                 "rootname":{'value':'root', 'option':''},
                 ## record 4a - output options
                 "iOptPrnt":{'value':0, 'option':''},
                 "iOptDigital":{'value':1, 'option':''},
                 "iOptExcelS":{'value':0, 'option':''},
                 "iOptExcelM":{'value':0, 'option':''},
                 "iOptRad":{'value':1, 'option':''},
                 "nwskip":{'value':1, 'option':''},
                 ## record 4b - model options
                 "iIOPmodel":{'value':2, 'option':''},
                 "iSkyRadModel":{'value':1, 'option':''},
                 "iSkyIrradModel":{'value':0, 'option':''},
                 "iChl":{'value':2, 'option':''},
                 "iCDOM":{'value':3, 'option':''},
                 ## record 5a - number of components
                 "ncomp":{'value':4, 'option':''},
                 "nconc":{'value':4, 'option':''},
                 ## record 6 - number of wavelengths
                 "nwave":{'value':0, 'option':''},
                 ## record 7 - inelastic scattering
                 "ibiolum":{'value':0, 'option':''},
                 "ichlfl":{'value':0, 'option':''},
                 "icdomfl":{'value':0, 'option':''},
                 "iraman":{'value':0, 'option':''},
                 "icompchl":{'value':2, 'option':''},
                 ## record 8a - sky model
                 "iflagsky":{'value':1, 'option':''},
                 "nsky":{'value':5, 'option':''},
                 "suntheta":{'value':30, 'option':''},
                 "sunphi":{'value':0.0, 'option':''},
                 "C":{'value':1.25, 'option':''},
                 "rsky":{'value':0.3, 'option':''},
                 "Edtotal":{'value':0.3, 'option':''},
                 "cloud":{'value':0., 'option':''},
                 "hour":{'value':0., 'option':''},
                 ## record 8b - atmospheric conditions
                 "jday":{'value':180., 'option':''},
                 "rlat":{'value':52., 'option':''},
                 "rlon":{'value':3., 'option':''},
                 "pres":{'value':29.92, 'option':''},
                 "am":{'value':1., 'option':''},
                 "rh":{'value':90., 'option':''},
                 "wv":{'value':2.5, 'option':''},
                 "vi":{'value':15., 'option':''},
                 "wsm":{'value':5., 'option':''},
                 "ro3":{'value':-99, 'option':''},
                 ## record 9 - surface information
                 "windspd":{'value':0, 'option':''},
                 "refr":{'value':1.34, 'option':''},
                 "temp":{'value':20, 'option':''},
                 "salinty":{'value':35, 'option':''},
                 ## record 10 - bottom reflectance
                 "ibotm":{'value':0, 'option':''},
                 "rflbot":{'value':0.2, 'option':''},
                 ## record 11 - output depths
                 "iop":{'value':0, 'option':''},
                 "nznom":{'value':1, 'option':''},
                 "zetanom":{'value':[0], 'option':''},
                 ## record 12 - data files
                'PureWaterDataFile':{'value':"../data/H2OabDefaults_SEAwater.txt", 'option':''},
                'nac9Files':{'value':1, 'option':''},
                'ac9DataFile':{'value':"dummyac9.txt", 'option':''},
                'Ac9FilteredDataFile':{'value':"dummyFilteredAc9.txt", 'option':''},
                'HydroScatDataFile':{'value':"dummyHscat.txt", 'option':''},
                'ChlzDataFile':{'value':"dummyComp.txt", 'option':''},
                'CDOMDataFile':{'value':"dummyComp.txt", 'option':''},
                'RbottomFile':{'value':"dummyR.bot", 'option':''},
                'TxtDataFile':{'value':["dummydata.txt"], 'option':''},
                'IrradDataFile':{'value':"DummyIrrad.txt", 'option':''},
                #'IrradDataFile':{'value':"../data/IOPS/Sky_Irrad_Example_Eddir_Eddif_1500.txt", 'option':''},
                'S0biolumFile':{'value':"../data/MyBiolumData.txt", 'option':''}
                }

        ## reset with user values
        for iv,v in enumerate(rdefs):
            if v in kwargs:
                rdefs[v]['value'] = kwargs[v]

        ## set wavelengths
        if 'wavel' in kwargs:
            wl = kwargs['wavel']
            if type(wl) != list: wl = [wl]
            if len(wl) == 1:
                rdefs['nwave']['value']=0
            else:
                rdefs['nwave']['value']=len(wl)-1
            rdefs['wavel']={'value':wl}
        else:
            waves = np.arange(min_wave-step,max_wave+step*2, step=step)
            nwaves = len(waves)-1
            rdefs['nwave']['value']=nwaves
            rdefs['wavel']={'value':list(waves)}

        ## add component specific values records 5b-h
        ## record 5b - component concentrations
        compvars = ['compconc{}'.format(i+1) for i in range(0,rdefs['nconc']['value'])]
        for ci, cv in enumerate(compvars):
            if ci == 0: compconc = []
            compconc.append('{}'.format(0 if not cv in kwargs else kwargs[cv]))
        rdefs['compconc']=','.join(compconc)

        ## record 5c - specific absorbtion parameters
        ## QV 2018-04-23 added itype here instead of with the phase functions
        defaults_abs = [{'itype':0,'iastropt':1, 'astarRef':-999, 'astar0':0.0,'asgamma':0.0},
                        {'itype':0,'iastropt':0, 'astarRef':-999, 'astar0':0.0,'asgamma':0.0},
                        {'itype':0,'iastropt':4, 'astarRef':440, 'astar0':0.06,'asgamma':0.014},
                        {'itype':0,'iastropt':0, 'astarRef':-999, 'astar0':0.0,'asgamma':0.0}]

        for i in range(0,rdefs['nconc']['value']):
            defaults = defaults_abs[i]
            for ip,par in enumerate(defaults):
                cpar = '{}{}'.format(par, i+1)
                cval='{}'.format(defaults[par] if not cpar in kwargs else kwargs[cpar])
                if ip == 0: curline = [] #"{}".format(i)
                curline.append('{}'.format(cval))
            rdefs['abspar{}'.format(i+1)]=', '.join(curline)

        ## record 5d - specific absorbtion file names
        defaults = {'astarfile':''}
        default_files = {'astarfile1':"../data/H2OabDefaults_SEAwater.txt",
                         'astarfile2':"../data/defaults/astarchl.txt",
                         'astarfile3':"dummyastar.txt",
                         'astarfile4':"../data/defaults/astarmin_average.txt"}
        for i in range(0,rdefs['nconc']['value']):
            for ip,par in enumerate(defaults):
                cpar = '{}{}'.format(par, i+1)
                cval='{}'.format(default_files[cpar] if not cpar in kwargs else kwargs[cpar])
                if ip == 0: curline = []
                curline.append('{}'.format(cval))
            rdefs['absfile{}'.format(i+1)]=', '.join(curline)

        ## record 5e - specific scattering parameters
        defaults_sca = [{'ibstropt':0, 'bstarRef':-999, 'bstar0':0.0,'coef1':0.0,'coef2':0.0,'coef3':0.0},
                        {'ibstropt':-1, 'bstarRef':-999, 'bstar0':0.0,'coef1':0.0,'coef2':0.0,'coef3':0.0},
                        {'ibstropt':-1, 'bstarRef':-999, 'bstar0':0.0,'coef1':0.0,'coef2':0.0,'coef3':0.0},
                        {'ibstropt':0, 'bstarRef':-999, 'bstar0':0.0,'coef1':0.0,'coef2':0.0,'coef3':0.0}]

        for i in range(0,rdefs['ncomp']['value']):
            defaults = defaults_sca[i]
            for ip,par in enumerate(defaults):
                cpar = '{}{}'.format(par, i+1)
                cval = defaults[par] if not cpar in kwargs else kwargs[cpar]
                if (par == "ibstropt") & (cval == 1):
                    defaults['bstarRef']=550
                    defaults['bstar0']=0.3
                    defaults['coef1']=1
                    defaults['coef2']=0.62

                if ip == 0: curline = []
                curline.append('{}'.format(cval))
            rdefs['scapar{}'.format(i+1)]=', '.join(curline)

        ## record 5f - specific scattering file names
        defaults = {'bstarfile':''}
        default_files = {'bstarfile1':"../data/H2OabDefaults_SEAwater.txt",
                         'bstarfile2':"../data/defaults/bstarchl.txt",
                         'bstarfile3':"bstarDummy.txt",
                         'bstarfile4':"../data/defaults/bstarmin_average.txt"}
        for i in range(0,rdefs['ncomp']['value']):
            for ip,par in enumerate(defaults):
                cpar = '{}{}'.format(par, i+1)
                cval='{}'.format(default_files[cpar] if not cpar in kwargs else kwargs[cpar])
                if ip == 0: curline = []
                curline.append('{}'.format(cval))
            rdefs['scafile{}'.format(i+1)]=', '.join(curline)

        ## record 5g - type and phase function of concentrations
        ## QV 2018-04-23 removed itype here
        defaults_pha = [{'ibbopt':0, 'bbfrac':0, 'BfrefPL':0, 'Bf0PL':0, 'BfmPL':0},
                        {'ibbopt':0, 'bbfrac':0, 'BfrefPL':0, 'Bf0PL':0, 'BfmPL':0},
                        {'ibbopt':0, 'bbfrac':0, 'BfrefPL':0, 'Bf0PL':0, 'BfmPL':0},
                        {'ibbopt':0, 'bbfrac':0, 'BfrefPL':0, 'Bf0PL':0, 'BfmPL':0}]

        for i in range(0,rdefs['ncomp']['value']):
            defaults = defaults_pha[i]
            for ip,par in enumerate(defaults):
                cpar = '{}{}'.format(par, i+1)
                cval='{}'.format(defaults[par] if not cpar in kwargs else kwargs[cpar])
                if ip == 0: curline = []
                #if par == 'itype': continue
                curline.append('{}'.format(cval))
            rdefs['phapar{}'.format(i+1)]=', '.join(curline)

        ## record 5h - phase function file names
        defaults = {'pfname':''}
        default_files = {'pfname1':"pureh2o.dpf",
                         'pfname2':"Case1Large.dpf",
                         'pfname3':"isotrop.dpf",
                         'pfname4':"avgpart.dpf"}
        for i in range(0,rdefs['ncomp']['value']):
            for ip,par in enumerate(defaults):
                cpar = '{}{}'.format(par, i+1)
                cval='{}'.format(default_files[cpar] if not cpar in kwargs else kwargs[cpar])
                if ip == 0: curline = []
                curline.append('{}'.format(cval))
            rdefs['phafile{}'.format(i+1)]=', '.join(curline)

        ## add wavelength records 6b
        if rdefs['nwave']['value']==0:
            defaults = {'wavel':550, 'areset':0.1, 'breset':0.25}
            for ip,par in enumerate(defaults):
                cval='{}'.format(defaults[par] if not par in kwargs else kwargs[par])
                if ip == 0: curline = []
                curline.append('{}'.format(cval))
            rdefs['waverec{}'.format(1)]=','.join(curline)
        elif rdefs['nwave']['value']>0:
            if 'wavel' in kwargs:
                cval = kwargs['wavel']
            else:
                cval = rdefs['wavel']['value']

            if type(cval) != list: cval = [cval]
            line, wline = '', 1
            for wi, wave in enumerate(cval):
                line += '{}, '.format(wave)
                if ((wi+1) % 12 == 0) | (wi == len(cval)-1):
                    rdefs['waverec{}'.format(wline)]=line
                    line = ''
                    wline += 1

        ## add sky model records 8
        if rdefs["iflagsky"]['value'] == 1:
            rdefs['nsky']['value']=5
            skytags = ["iflagsky",'nsky',"suntheta",'sunphi','C','rsky','Edtotal']
        elif rdefs["iflagsky"]['value'] == 2:
            rdefs['nsky']['value']=3
            skytags = ["iflagsky",'nsky',"suntheta",'sunphi','cloud']
        elif rdefs["iflagsky"]['value'] == 3:
            rdefs['nsky']['value']=3
            skytags = ["iflagsky",'nsky',"hour",'cloud','sunphi']

        cval=[rdefs[cv]['value'] for cv in skytags]
        rdefs['skyrec{}'.format(1)]=', '.join([str(c) for c in cval])

        ## add depth records 11
        rdefs['nznom']['value']=len(rdefs['zetanom']['value'])
        cval = [rdefs[r]['value'] for r in ['iop','nznom']]
        zout=1
        for zi, z in enumerate(rdefs['zetanom']['value']):
            cval.append(z)
            if (len(cval) == 12) | (zi == len(rdefs['zetanom']['value'])-1):
                rdefs['depthrec{}'.format(zout)]=', '.join([str(c) for c in cval])
                rdefs['depthrec{}'.format(zout)]+=','
                zout+=1
                cval = []
        #cval += rdefs['zetanom']['value']

        ## add Data files record 12
        datafiles = [""]*rdefs['ncomp']['value']
        for i in range(0,rdefs['ncomp']['value']):
            cv = 'TxtDataFile{}'.format(i+1)
            datafiles[i] = "dummydata.txt" if cv not in kwargs else kwargs[cv]
        rdefs['TxtDataFile']['value'] = '\n'.join(datafiles)

        ## add custom irrad
        for k in ['IrradDataFile']:
            if k in kwargs:
                rdefs[k]['value'] = kwargs[k]

        return(rdefs)
