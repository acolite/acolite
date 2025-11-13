## plane_fit
## function to perform a plane fit to a dataset along given x and y dimensions
## based on code from Arthur Coqué used to resample the 5000 m MSI angles grids
## changed here to get x and y as inputs, by default already as mesh
##
## written by Quinten Vanhellemont, RBINS
## 2025-11-13
## modifications:

def plane_fit(geom_x, geom_y, geom_grid, mesh = True):
    import numpy as np
    import scipy

    if not mesh:
        ii, jj = np.meshgrid(geom_x, geom_y, indexing='ij')
        arr = np.vstack([ii.ravel(), jj.ravel(), geom_grid.ravel()])
    else:
        arr = np.vstack([geom_x.ravel(), geom_y.ravel(), geom_grid.ravel()])

    points = arr[:, ~np.isnan(arr).any(axis=0)]  # (3, N)
    centroid = np.mean(points, axis=1)
    A = points - centroid.reshape(-1, 1)
    U, *_ = scipy.linalg.svd(A, full_matrices=False, overwrite_a=True, check_finite=False)
    normal = U[:, 2]
    d = normal[0] * centroid[0] + normal[1] * centroid[1] + normal[2] * centroid[2]
    return(normal, d)
