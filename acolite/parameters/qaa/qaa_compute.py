## function to compute Lee's QAA with DVDZ merging
## QV 2021-02-15

def qaa_compute(qaa_in, sza = 0, satellite = None,
                qaa_wave = [443, 490, 560, 665],qaa_coef = None,
                compute_zeu_lee = False, k_atag = 'v6_a_490',k_bbptag = 'v6_bbp_490'):

    import time
    import numpy as np
    import acolite as ac

    if qaa_coef is None: qaa_coef = ac.parameters.qaa.qaa_coef()
    t0 = time.time()

    ## indices for v5 or v6 cal
    v5_idx = np.where(qaa_in[665] < 0.0015)
    v6_idx = np.where(qaa_in[665] >= 0.0015)

    ## run through bands
    qaa_data = {}
    for iw, k in enumerate(qaa_in):
        ## do band shifting with coefficients provided by Dimitry
        if qaa_coef['spectral_shift']:
            ctag = 'Rrs{}_{}_a'.format(k, satellite)
            qaa_in[k] = (qaa_coef[ctag][0]*qaa_in[k]*qaa_in[k]) + (qaa_coef[ctag][1]*qaa_in[k]) + (qaa_coef[ctag][2])

        ## step 0
        ## compute Lee's subsurface rrs
        rname = 'rrs_{}'.format(k)
        qaa_data[rname] = qaa_in[k]/(0.52+(1.7*qaa_in[k]))

        ## step 1
        ## compute Lee's u
        uname = 'u_{}'.format(k)
        qaa_data[uname] = (-1.0*qaa_coef['g'][0] + np.sqrt((np.power(qaa_coef['g'][0],2.)+(4.*qaa_coef['g'][1]*qaa_data[rname]))))/(2.*qaa_coef['g'][1])

    ## step 2
    ## get rhow and a560 according to QAA5
    rhow_v5=np.log10((qaa_data['rrs_443']+qaa_data['rrs_490'])/\
                        (qaa_data['rrs_560']+(5.*qaa_data['rrs_665']*(qaa_data['rrs_665']/qaa_data['rrs_490']))))
    qaa_data['v5_a_560'] = qaa_coef['aw'][2]+np.power(10.,(qaa_coef['h'][2]+(qaa_coef['h'][1]*rhow_v5)+qaa_coef['h'][0]*np.power(rhow_v5,2.)))

    ## get rhow and a665 according to QAA6
    rhow_v6=qaa_data['rrs_665']/(qaa_data['rrs_443']+qaa_data['rrs_490'])
    qaa_data['v6_a_665']= qaa_coef['aw'][3]+qaa_coef['k'][0]*np.power(rhow_v6,qaa_coef['k'][1])

    ## step 3
    qaa_data['v5_bbp_560'] = ((qaa_data['u_560'] * qaa_data['v5_a_560']) / (1.0 - qaa_data['u_560'])) - qaa_coef['bbw'][2]
    qaa_data['v6_bbp_665'] = ((qaa_data['u_665'] * qaa_data['v6_a_665']) / (1.0 - qaa_data['u_665'])) - qaa_coef['bbw'][3]

    ## step 4 QAAv6
    ratio = qaa_data['rrs_443']/qaa_data['rrs_560']
    Y = qaa_coef['l'][0] * (1.0 - qaa_coef['l'][1] * np.exp(qaa_coef['l'][2]*ratio))

    ## compute a QAA5 & QAA6
    for version in ['5','6']:
        ## reference wavelengths
        refw = {'5':555, '6':670}[version]
        refn = {'5':560, '6':665}[version]
        for iw,wave in enumerate(qaa_wave):
            ## get dataset tags
            bbptag = 'v{}_bbp_{}'.format(version,wave)
            atag = 'v{}_a_{}'.format(version,wave)
            utag = 'u_{}'.format(wave)
            kdtag = 'v{}_Kd_{}'.format(version,wave)

            ## compute bbp and q
            if not (((version == '5') & (wave == 560)) or ((version == '6') & (wave == 665))):
                ##step 5 - get bbp
                qaa_data[bbptag] = qaa_data['v{}_bbp_{}'.format(version,refn)] * \
                                       np.power(float(refw)/float(wave),Y)
                ##step 6 - get a
                qaa_data[atag] = ((1.0 - qaa_data[utag]) * (qaa_coef['bbw'][iw]+qaa_data[bbptag])) / qaa_data[utag]

            ## compute Kd
            #qaa_data[kdtag] = qaa_kd(qaa_data[atag], qaa_data[bbptag], qaa_coef)
            v=qaa_coef['m'][1]*(1.-qaa_coef['m'][2]*np.exp(-qaa_coef['m'][3]*qaa_data[atag]))
            qaa_data[kdtag]=qaa_coef['m'][0]*qaa_data[atag]+v*qaa_data[bbptag]

    ## get switched datasets
    for iw,wave in enumerate(qaa_wave):
        for qaapar in ['a','bbp','Kd']:
            t5 = 'v5_{}_{}'.format(qaapar,wave)
            t6 = 'v6_{}_{}'.format(qaapar,wave)
            #tsw = '{}{}_vw'.format(qaapar,wave)
            #qaa_data[tsw] = qaa_data[t6]
            #if len(v5_idx[0]>0): qaa_data[tsw][v5_idx]=qaa_data[t5][v5_idx]
            if len(v5_idx[0]>0): qaa_data[t6][v5_idx]=qaa_data[t5][v5_idx]

    ## parameters for KdPAR 0-1m Nechad KdPARv2
    qaa_data['v5_KdPAR_Nechad']=1.2529*np.log(1.+qaa_data['v5_Kd_490'])+0.1127
    qaa_data['v6_KdPAR_Nechad']=1.2529*np.log(1.+qaa_data['v6_Kd_490'])+0.1127

    ## get Zeu
    c1=[-0.057,0.482,4.221]
    c2=[0.183,0.702,-2.567]
    alpha=[0.090,1.465,-0.667]

    stheta = np.sin(sza*(np.pi/180))
    ctheta = np.cos(sza*(np.pi/180))
    tau = -np.log(0.01) # Zeu

    k1 = (c1[0] + c1[1]*np.sqrt(qaa_data[k_atag]) + c1[2]*qaa_data[k_bbptag])*(1+alpha[0]*stheta);
    k2 = (c2[0] + c2[1]*       (qaa_data[k_atag]) + c2[2]*qaa_data[k_bbptag])*(alpha[1]+alpha[2]*ctheta)

    # Kpar, Lee et al. 2007
    z = 1.
    qaa_data['v6_KPAR_Lee'] = k1+(k2)/np.sqrt(z+1.)
    qaa_data['v6_Zeu_Lee'] = 4.6/qaa_data['v6_KPAR_Lee']

    ## slow (root solving)
    if compute_zeu_lee:
        y1 = (k1*k1 - k2*k2 - 2.0*tau*k1)/(k1*k1)
        y2 = (tau*tau - 2.0*tau*k1)/(k1*k1)
        y3 = (tau*tau)/(k1*k1)
        zeu = qaa_data[atag]*0.0
        zeu[:]=np.nan

        val = np.where(np.isfinite(y1) & np.isfinite(y2) & np.isfinite(y3))
        for i,i1 in enumerate(val[0]):
            i2=val[1][i]
            zeu[i1,i2] = np.real(np.roots((y1[i1,i2],y2[i1,i2],y3[i1,i2])))[1]
        zeu[zeu<0]=np.nan
        qaa_data['v6_Zeu_Lee_roots'] = zeu

    t1 = time.time()
    print('Finished QAA computations in {:.1f} seconds'.format(t1-t0))

    return(qaa_data)
