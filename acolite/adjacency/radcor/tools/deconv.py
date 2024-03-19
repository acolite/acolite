## def deconv
## perform deconvolution
##
## written by Alexandre Castagna, UGent (R version)
##            Quinten Vanhellemont, RBINS (Python version)
##
## for the RAdCor project 2023-2024
##
## modifications: 2024-03-12 (QV) added as function

def deconv(x_f, otf, id_psf, x_f_dim, keep_edges = False):
    import numpy as np
    res = np.fft.ifftn(x_f / otf).real
    if not keep_edges:
        idr = np.array(range(x_f_dim[0] - id_psf, x_f_dim[0]))
        idr = np.append(idr, range(0, x_f_dim[0] - id_psf))
        idc = np.array(range(x_f_dim[1] - id_psf, x_f_dim[1]))
        idc = np.append(idc, range(0, x_f_dim[1] - id_psf))
    else:
        idr = np.array(range(0, x_f_dim[0] - id_psf * 2))
        idc = np.array(range(0, x_f_dim[0] - id_psf * 2))
    return(res[idr][:, idc])
