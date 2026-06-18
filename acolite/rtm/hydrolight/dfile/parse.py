## hydrolight.dfile.parse
## reads Hydrolight dfile and formats to quads
## computes Rrs
##
## written by Quinten Vanhellemont, RBINS
## 2018-04-24
## modifications: 2021-09-27 (QV) reversed order of nphi and nmu for reshape
##                2021-09-28 (QV) track azi and view, transpose 2d data and add 360 azimuth
##                                new reshape for depth resolved 3d data, added 360 azimuth
##                2021-10-05 (QV) fill 0 degree cap with same values
##                2026-06-16 (QV) new function, renamed wavelengths to wavelength, copy wavelength as array
##                                renamed from dfile_import
def parse(inp):
    import acolite as ac
    import numpy as np

    if type(inp) is dict:
        hd = {k: inp[k] for k in inp}
    else:
        ## read data into dict
        hd = ac.rtm.hydrolight.dfile.read(inp)

    from numpy import asarray
    #empty = ndarray((hd['nmu'],hd['nphi'],hd['nwave']))

    data = {'wavelength': np.asarray(hd['wavelength']),
            'phi': np.asarray(hd[hd['wavelength'][0]]['phi']['data']),
            'fmu': np.asarray(hd[hd['wavelength'][0]]['fmu']['data'])}
    data['view'] = np.flip(data['fmu'] * (180/np.pi))
    data['azi'] = np.append(data['phi'], data['phi'][-1] + np.diff(data['phi'])[-1])

    for wi, wave in enumerate(hd['wavelength']):
        if wave not in data: data[wave] = {}

        for par in ['Ed',
                    'Lw_air','Lrefl_air', 'Lu_air','Lsky_air',
                    'Lu_water', 'Ld_water']:
            simpar = None
            reshape = None #(hd['nmu'],hd['nphi'])

            ## 1d depth
            if par in ['Ed']:
                reshape = (hd['nz'])
            if par == 'Ed': simpar = 'Ed'

            ## 2d quads
            if par in ['Lw_air','Lrefl_air', 'Lu_air', 'Lsky_air']:
                #reshape = (hd['nmu'],hd['nphi'])
                reshape = (hd['nphi'],hd['nmu'])

            if par == 'Lw_air': simpar = 'RADMa'
            if par == 'Lrefl_air': simpar = 'RAD0Ma'
            if par == 'Lu_air': simpar = ['RADMa','RAD0Ma']
            if par == 'Lsky_air': simpar = 'radsky'

            ## 3d quads with depth
            if par in ['Lu_water', 'Ld_water']:
                #reshape = (hd['nmu'],hd['nphi'],hd['nz'])
                #reshape = (hd['nphi'],hd['nmu'],hd['nz'])
                reshape = hd['nz'],hd['nphi'],hd['nmu']

            ## diffuse up and down radiances in water
            if par == 'Lu_water': simpar = 'RADMz'
            if par == 'Ld_water': simpar = ['RADPz', 'RAD0Pz']

            if simpar == None: continue

            ## get data from hd dict
            if type(simpar) == list:
                for si, sim in enumerate(simpar):
                    if si == 0:
                        cur = asarray(hd[wave][sim]['data'])
                    else:
                        cur += asarray(hd[wave][sim]['data'])
            else:
                cur = asarray(hd[wave][simpar]['data'])

            if reshape == None:
                data[wave][par] = cur
            elif type(reshape)==int:
                data[wave][par] = cur
            elif len(reshape)==2:
                #data[wave][par] = cur.reshape(reshape[0],reshape[1])

                ## repeat element for 360 azi
                pd = cur.reshape(reshape[0],reshape[1]).T
                pd = np.hstack((pd, np.atleast_2d(pd[:, 0]).T))

                ## replace zeros in nadir view
                pd[hd['nmu']-1, :] = pd[hd['nmu']-1, 0]

                data[wave][par] = pd

            elif len(reshape)==3:
                #data[wave][par] = cur.reshape(reshape[0],reshape[1], reshape[2])
                pdi = cur.reshape(reshape[0],reshape[1], reshape[2])
                for i in range(reshape[0]):
                    ## repeat element for 360 azi
                    pd = pdi[i,:,:].T
                    pd = np.hstack((pd, np.atleast_2d(pd[:, 0]).T))

                    ## replace zeros in nadir view
                    pd[hd['nmu']-1, :] = pd[hd['nmu']-1, 0]

                    if i == 0:
                        pdo = pd
                    else:
                        pdo = np.dstack((pdo, pd))
                data[wave][par] = pdo
            else:
                data[wave][par] = cur

        data[wave]['Rrs']=data[wave]['Lw_air']/data[wave]['Ed'][0]
    return(data)
