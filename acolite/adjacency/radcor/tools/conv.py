## def conv
## perform convolution
##
## written by Alexandre Castagna, UGent (R version)
##            Quinten Vanhellemont, RBINS (Python version)
##
## for the RAdCor project 2023-2024
##
## modifications: 2024-03-12 (QV) added as function

def conv(x_f, otf, id_psf, x_f_dim, keep_edges = False):
    import numpy as np
    res = np.fft.ifftn(otf * x_f).real
    if not keep_edges:
        idr = np.array(range(id_psf, x_f_dim[0]))
        idr = np.append(idr, range(0, id_psf))
        idc = np.array(range(id_psf, x_f_dim[1]))
        idc = np.append(idc, range(0, id_psf))
    else:
        idr = np.array(range(id_psf * 2, x_f_dim[0]))
        idc = np.array(range(id_psf * 2, x_f_dim[1]))
    return(res[idr][:, idc])
