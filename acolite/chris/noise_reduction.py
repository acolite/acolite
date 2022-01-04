## do CHRIS noise reduction/destriping
## original Python Code by Héloïse Lavigne
## function for ACOLITE processing QV 2021-06-09
## modifications: 2021-12-31 (QV) skip TOA radiances when creating the RTOA dataset
##                2022-01-04 (QV) added netcdf compression

def noise_reduction(ncf, rename=True,
                         netcdf_compression=False,
                         netcdf_compression_level=4,
                         netcdf_compression_least_significant_digit=None):
    import numpy as np
    from scipy.ndimage import gaussian_filter
    import acolite as ac

    ## function from spectral
    def spectral_angles(data, members):
        '''Calculates spectral angles with respect to given set of spectra.
        Arguments:
            `data` (:class:`numpy.ndarray` or :class:`spectral.Image`):
                An `MxNxB` image for which spectral angles will be calculated.
            `members` (:class:`numpy.ndarray`):
                `CxB` array of spectral endmembers.
        Returns:
            `MxNxC` array of spectral angles.
        Calculates the spectral angles between each vector in data and each of the
        endmembers.  The output of this function (angles) can be used to classify
        the data by minimum spectral angle by calling argmin(angles).
        '''
        assert members.shape[1] == data.shape[2], \
            'Matrix dimensions are not aligned.'

        m = np.array(members, np.float64)
        m /= np.sqrt(np.einsum('ij,ij->i', m, m))[:, np.newaxis]

        norms = np.sqrt(np.einsum('ijk,ijk->ij', data, data))
        dots = np.einsum('ijk,mk->ijm', data, m)
        dots = np.clip(dots / norms[:, :, np.newaxis], -1, 1)
        return np.arccos(dots)

    ## get datasets and attributes
    datasets = ac.shared.nc_datasets(ncf)
    gatts = ac.shared.nc_gatts(ncf)

    ## read data
    ds_att = {}
    RTOA = None
    for di, ds in enumerate(datasets):
        if 'rhot_' not in ds: continue
        tmp, att = ac.shared.nc_data(ncf, ds, attributes=True)
        ds_att[ds] = att
        if RTOA is None:
            RTOA = tmp
        else:
            RTOA = np.dstack((RTOA, tmp))

    ## do some masking
    RTOA = np.where(np.isnan(RTOA), -5, RTOA)
    RTOA = np.where(RTOA <= 0.01, np.nan, RTOA)

    ## get array shape
    nrow, ncol, nwl = RTOA.shape

    ## interband calibration already done
    RTOAcal1 = RTOA * 1.0

    ### 2. apply a column smoothing.
    ### 2.A Dropout correction (Gomez-Chova et al. 2008 Applied Optics)
    Dall = np.ones([nrow, ncol, nwl])*np.nan
    Deven = np.ones([nrow, ncol, nwl])*np.nan

    for w in range(nwl):
        for p in range(ncol-1):
            Dall[:,p,w] = (RTOAcal1[:,p,w] - RTOAcal1[:,p+1,w])**2

    for w in range(nwl):
        for p in np.arange(1, ncol-2, 2):
            Deven[:,p,w] = (RTOAcal1[:,p,w] - RTOAcal1[:,p+2,w])**2

    DROPOUT = np.ones([nrow, ncol, nwl])*0

    for w in range(nwl):
        for b in np.arange(0, nrow, 1):
            delta = Dall[b,:,w]/Deven[b,:,w]
            mdelta = np.nanmedian(delta)
            if mdelta >= 1.5:
                DROPOUT[b,np.arange(0, ncol, 2) ,w] = 1

    ### Correct dropout pixels
    DD = np.where(DROPOUT == 1)
    ndrop = len(DD[0])

    RTOAcal1d = RTOAcal1 * 1.0
    for d in range(ndrop):
        irow=DD[0][d]
        icol=DD[1][d]
        iwl=DD[2][d]
        if np.isnan(RTOAcal1[irow, icol, iwl]) == False:
            wmin = max(iwl-2, 0)
            wmax = min(iwl+2, 61)
            Wup = 0
            Wdown = 0
            for w in np.arange(wmin, wmax, 1):
                if irow != 0 and w != iwl:
                    Wup = Wup + (RTOAcal1[irow, icol, iwl] - RTOAcal1[irow-1, icol, w])**2
                if irow != nrow-1 and w != iwl:
                    Wdown = Wdown + (RTOAcal1[irow, icol, iwl] - RTOAcal1[irow+1, icol, w])**2
            Wup = Wup**(-1/2) if Wup > 0 else 0
            Wdown = Wdown**(-1/2) if Wdown > 0 else 0
            Wupc = Wup / (Wup + Wdown) if Wup+Wdown != 0 else 0.5
            Wdownc = Wdown / (Wup + Wdown) if Wup+Wdown != 0 else 0.5
            if irow != 0 and irow != nrow-1:
                RTOAcal1d[irow, icol, iwl] = RTOAcal1[irow+1, icol, iwl]*Wupc + RTOAcal1[irow-1, icol, iwl]*Wdownc
            if irow == 0 :
                RTOAcal1d[irow, icol, iwl] = RTOAcal1[irow+1, icol, iwl]*Wdownc
            if irow == nrow-1:
                RTOAcal1d[irow, icol, iwl] = RTOAcal1[irow-1, icol, iwl]*Wupc

    ### 2.B Vertical striping correction (Gomez-Chova et al. 2008 Applied Optics)
    ### identification of edges in the image (comparison of spectrum)
    DISTmat = np.ones([nrow, ncol])
    for i in range(ncol-1):  ##column
        for j in range(nrow):  ##ligne
            a = np.ones([1,1,62])
            a[0,0,:] = RTOAcal1d[j, i, :]
            b = np.ones([1,62])
            b[0,:] = RTOAcal1d[j,i+1,:]
            R = spectral_angles(a, b)
            DISTmat[j,i] = R[0][0][0]


    DISTmat[:,ncol-1] = DISTmat[:,ncol-2]

    DISTvect = np.nanquantile(DISTmat, 0.8, axis=1)


    ## si un pixel est identifié comme "high difference on masque toute la ligne)
    LIMIT = np.nanquantile(DISTvect, 0.6)
    MASK = np.ones([nrow, ncol])
    for i in range(nrow):
        if DISTvect[i] >= LIMIT and np.isnan(DISTvect[i])==False:
            MASK[i,:] = np.nan
        if np.isnan(DISTvect[i]):
            MASK[i,:] = np.nan

    ## verifie si le mask n'est ne contient pas uniquement des NANs, dans ce cas le remettre automatiquement à 1
    NNA = len(np.where(np.isnan(MASK))[1])
    if NNA >= 10000 :
        MASK = np.ones([nrow, ncol])


    ### apply vertical striping
    A = np.ones([nwl, ncol])*np.nan
    for w in range(nwl):
        TMP = RTOAcal1d[:,:,w]*MASK
        A[w,:] = np.nanmean(TMP,axis=0)
    B = np.log10(A)
    C = np.ones([nwl, ncol])*np.nan
    for w in range(nwl):
        C[w,:] = gaussian_filter(B[w,:], sigma=2)
    D=C-B
    E = 10**D

    ## final array
    RTOAcal2 = np.ones([nrow, ncol, nwl])*np.nan
    for w in range(nwl):
        for i in range(ncol):
            RTOAcal2[:,i,w] = RTOAcal1d[:,i,w]*E[w,i]

    ## output result
    if rename:
        ofile = ncf.replace('_L1R.nc', '_L1R_NR.nc')
        new = True
    else:
        ofile = '{}'.format(ncf)
        new = False
    dix = 0
    for di, ds in enumerate(datasets):
        if 'rhot_' in ds:
            ac.output.nc_write(ofile, ds, RTOAcal2[:,:,dix], dataset_attributes=ds_att[ds], attributes=gatts, new=new,
                               netcdf_compression=netcdf_compression, netcdf_compression_level=netcdf_compression_level,
                               netcdf_compression_least_significant_digit=netcdf_compression_least_significant_digit)
            dix += 1
        else:
            tmp, att = ac.shared.nc_data(ncf, ds, attributes=True)
            ac.output.nc_write(ofile, ds, tmp, dataset_attributes=att, attributes=gatts, new=new,
                               netcdf_compression=netcdf_compression, netcdf_compression_level=netcdf_compression_level,
                               netcdf_compression_least_significant_digit=netcdf_compression_least_significant_digit)
        new = False
    return(ofile)
