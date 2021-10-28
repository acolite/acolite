## make the atmospheric weigths
## adapted from AC R code v4
## QV 2020-01-29
## 2021-07-06 (QV) updated with new cumulative distribution prediction and ex/resolution now in km

def w_kernel(coef_aer, coef_ray, tray, taer, ex = 0, res = 0.03, pressure=1013.25):
    import acolite as ac
    import numpy as np
    from scipy.misc import derivative

    # Construct the the weighted cumulative function
    def fun_grad(r, pressure=1013.25):
        fa = ac.adjacency.acstar3.pred_annular_cdf(r, coef_aer, pressure=1) # aerosol cdf only at normal pressure
        fr = ac.adjacency.acstar3.pred_annular_cdf(r, coef_ray, pressure=pressure/1013.25) # normalize pressure to pred_annular_cdf
        ft = (tray * fr + taer * fa) / (tray + taer)
        if type(ft) is np.ndarray:
            ft[np.isnan(ft)] = 0
        else:
            if np.isnan(ft): ft = 0
        return(ft)

    ## if ex == 0 optimize for 60% of diffuse transmittance
    if ex == 0:
        import scipy.optimize
        ## function to optimize
        def of(x, percentage = 0.6):
            if type(x) != np.ndarray:
                x = np.asarray([x])
            return(abs(fun_grad(x) - percentage))
        ex = scipy.optimize.brent(of)[0]
        print('Optimized ex: {:.2f}km'.format(ex))

    if ex < (res + res / 2):
        print("Radial extent has to be equal of larger than res+res/2")

    ## some np manipulation to get the xy grid
    ## there are probably better functions - but hard to find with no internet :-D
    xvec = np.arange(res/2, ex, res)
    xvec = np.hstack((-1*np.flip(xvec), xvec))
    xvec = xvec.reshape(len(xvec),1)
    xvec = np.tile(xvec, (1, len(xvec)))

    ## in fact xvec is yvec and xvec needs to be turned
    yvec = xvec * 1
    xvec = np.rot90(xvec)

    ## don't flatten now - we can keep 2D
    if False:
        yvec = yvec.flatten()
        xvec = xvec.flatten()

    ## in km
    rm = np.sqrt(np.power(xvec,2)+np.power(yvec,2))


    ##
    wm = derivative(fun_grad, rm, dx=0.001) / 2 / np.pi / rm
    wm = wm[1:,] + wm[:-1,]
    wm = wm[:,1:] + wm[:,:-1]
    wm = (res)**2 * wm / 4

    ## recompute center point weight
    rgt = np.sqrt((res)**2 / np.pi)
    xx,yy = np.where(wm == np.nanmax(wm))
    wm[xx, yy] = fun_grad(np.asarray([rgt]))

    return(wm, ex)
