## def nc_dim
## get dataset dimensions
##
## written by Alexandre Castagna, UGent (R version)
##            Quinten Vanhellemont, RBINS (Python version)
##
## for the RAdCor project 2023-2024
##
## modifications: 2024-03-12 (QV) added as function

def nc_dim(file):
    from netCDF4 import Dataset
    with Dataset(file) as nc:
        dim = len(nc.dimensions['y']), len(nc.dimensions['x'])
    return(dim)
