## move psf to corners and apply fft
## QV 2021-07-07 from AC code v0.9

def psf_otf(dim, psf, edge=True, fft=True):
    import numpy as np
    pdim = psf.shape

    idp = int((pdim[0] - 1) / 2)
    if edge: dim = dim[0]+idp * 2, dim[1]+idp * 2

    otf = np.zeros(dim)#+np.nan

    ## top left
    otf[0:idp+1, 0:idp+1] = psf[idp:, idp:]

    ## bottom left
    otf[(dim[0] - idp):, 0:idp+1] =  psf[:idp, idp:]

    ## bottom right
    otf[(dim[0] - idp):, (dim[1] - idp):] =  psf[:idp, :idp]

    ## top right
    otf[0:idp+1, (dim[1] - idp):] =  psf[idp:, :idp]

    ## apply 2d fft
    if fft:
        otf = np.fft.fftn(otf)

    return(otf)
