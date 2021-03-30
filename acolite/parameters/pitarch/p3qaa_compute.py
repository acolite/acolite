## function to compute Jaime's 3 band QAA
## QV 2021-02-15

def p3qaa_compute(sensor, b, g, r, cfg=None):
    import acolite as ac
    import numpy as np

    if cfg is None:
        cfg = ac.parameters.pitarch.p3qaa_coef()
        if sensor not in cfg:
            print('{} not in configuration'.format(sensor))
            return

    ## constants
    gamma=0.265
    m1=4.259
    m2=0.52
    m3=10.8
    ## QAA Rrs to u coefficients (Lee et al., (2002))
    g0=0.089
    g1=0.1245

    B = np.atleast_2d(b)
    G = np.atleast_2d(g)
    R = np.atleast_2d(r)

    Rrs = np.dstack((B,G,R))

#    ## R correction if off reasonable limits # 20210215 version
#    R_max=20*Rrs[:,:,1]**1.5
#    R_min=0.9*Rrs[:,:,1]**1.7
#    sub = np.where((Rrs[:,:,2]<R_min) | (Rrs[:,:,2]>R_max))
#    if len(sub[0])>0:
#        Rrs[sub[0],sub[1],2]=1.27*Rrs[sub[0],sub[1],1]**1.47+0.00018*(Rrs[sub[0],sub[1],2]/Rrs[sub[0],sub[1],1])**3.19
#    ## End of Rrs(670) correction if off reasonable limits

    ## R correction if off reasonable limits # 20210215 version
    R_max=20*Rrs[:,:,1]**1.5
    sub = np.where((Rrs[:,:,2]>R_max))
    if len(sub[0])>0:
        Rrs[sub[0],sub[1],2]=1.27*Rrs[sub[0],sub[1],1]**1.47+0.00018*(Rrs[sub[0],sub[1],2]/Rrs[sub[0],sub[1],1])**3.19
    sub = np.where((Rrs[:,:,2]<0))
    if len(sub[0])>0:
        ## I found this relationship by linear regression of Rrs(670) vs. Rrs(555). Nevertheless, this value is pretty irrelevant when the R band is small.
        Rrs[sub[0],sub[1],2]=0.182*Rrs[sub[0],sub[1],1]
    ## End of Rrs(670) correction if off reasonable limits


    ## Raman scattering correction
    bg_ratio = np.asarray([cfg[sensor]['bg_ratio'][k] for k in ['p4','p3','p2','p1','p0']])
    pscaled = np.poly1d(bg_ratio)
    al = np.asarray([cfg[sensor]['alpha'][k] for k in ['B', 'G', 'R']])
    b1 = np.asarray([cfg[sensor]['beta1'][k] for k in ['B', 'G', 'R']])
    b2 = np.asarray([cfg[sensor]['beta2'][k] for k in ['B', 'G', 'R']])
    ps = pscaled(B/G)
    for bi in range(3):
        RF = al[bi] * ps + b1[bi] * G ** b2[bi]
        Rrs[:,:,bi] = np.multiply(Rrs[:,:,bi], 1/(1+RF))
    rrs=Rrs/(0.52+1.7*Rrs) # rrs from Rrs
    ## End of Raman scattering correction

    ## Estimation of absorption at a reference band
    p_chi = np.asarray([cfg[sensor]['chi'][k] for k in ['p3','p2','p1','p0']])
    cscaled = np.poly1d(p_chi)
    aw = np.asarray([cfg[sensor]['aw'][k] for k in ['B', 'G', 'R']])
    wave = np.asarray([cfg[sensor]['center_wl'][k] for k in ['B', 'G', 'R']])
    i_l0=1 # 2 in Matlab
    l0=wave[i_l0]
    chi=np.log10(2*B/(G+5*R**2/B))
    al0=aw[i_l0]+10**cscaled(chi)
    ##End of estimation of absorption at a reference band

    bbw = np.asarray([cfg[sensor]['bbw'][k] for k in ['B', 'G', 'R']])
    u=((-g0+(g0**2+4*g1*rrs)**0.5)/(2*g1))
    bbpl0 = u[:,:,i_l0]*al0/(1-u[:,:,i_l0])-bbw[i_l0]
    eta=2*(1-1.2*np.exp(-0.9*pscaled(B/G)))
    bbp=Rrs*0
    for bi in range(3):
        bbp[:,:,bi] = bbpl0[:,:]*(l0/wave[bi])**eta[:,:]
    bb=bbp+bbw

    ## compute absorption
    a=(1-u)*bb/u
    ## If absorption is lower than that of water, it is set to water
    for bi in range(3):
        a[:,:,bi][a[:,:,bi]<aw[bi]] = aw[bi]
    ## Recalculation of bb to ensure closure
    bb = u*a/(1-u)

    ## compute Kd model (Lee et al., (2013))
    bbw_frac=bbw/bb
    Kd=a+(1-gamma*bbw_frac)*m1*(1-m2*np.exp(-m3*a))*bb #Spectral diffuse attenuation coefficient

    ## compute Secchi depth z_SD model (Lee et al., (2015))
    Kdidx = np.argsort(Kd, axis=2)[:,:,0]
    #z_SD_biased=1/(2.5*Kd)*np.log(np.abs(0.14-Rrs)/0.013)# % Lee et al. (2015)
    #return(z_SD_biased, Kdidx, Kd)

    #z_SD_biased = z_SD_biased[Kdidx[:,:,0].flatten()]
    z_SD_biased = np.zeros(R.shape)
    for bi in range(3):
        z_SD_biased[Kdidx==bi] = 1/(2.5*Kd[:,:,bi][Kdidx==bi])*np.log(np.abs(0.14-Rrs[:,:,bi][Kdidx==bi])/0.013)# % Lee et al. (2015)


    p_zsd = np.asarray([cfg[sensor]['coef_z_sd'][k] for k in ['p3','p2','p1','p0']])
    wscaled = np.poly1d(p_zsd)

    z_SD = wscaled(z_SD_biased)

    return({'a':a,'bb':bb,'Kd':Kd,'zSD':z_SD,'zSD_biased':z_SD_biased, 'eta':eta})
