## def extend
## select option to extend image when the edges need to be kept
##
## written by Alexandre Castagna, UGent (R version)
##            Quinten Vanhellemont, RBINS (Python version)
##
## for the RAdCor project 2023-2024
##
## modifications: 2024-03-12 (QV) added as function, compute average here

def extend(x_a, id_psf, x_a_dim, x_f_dim, method = 'average'):
    import acolite as ac
    import numpy as np

    ## Coordinates of original image in extended scene
    idr_o = (id_psf), (x_a_dim[0] + id_psf)
    idc_o = (id_psf), (x_a_dim[1] + id_psf)
    if (method == "average"):
        x_f = np.zeros(x_f_dim) + np.nanmean(x_a)
        x_f[idr_o[0]:idr_o[1], idc_o[0]:idc_o[1]] = x_a
    elif (method == "mirror"):
        x_f = ac.adjacency.radcor.extend.mirror(x_a, id_psf)
    elif (method == "nearest"):
        x_f = np.zeros(x_f_dim) + np.nan
        x_f[idr_o[0]:idr_o[1], idc_o[0]:idc_o[1]] = x_a
        x_f = ac.shared.fillnan(x_f)
    else:
        print('Error: Method {} not configured.'.format(method))
        return
    return(x_f)
