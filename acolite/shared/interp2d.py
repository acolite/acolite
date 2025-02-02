## class interp2d
## linear interpolator which returns NN for values outside the edge
## written by Quinten Vanhellemont, RBINS
## 2025-02-02
## modifications:
##

import numpy as np
import scipy.interpolate

class interp2d:
    def __init__(self, xy, v, outside_nn = True):
        self.lni = scipy.interpolate.LinearNDInterpolator(xy, v)
        self.outside_nn = outside_nn
        if self.outside_nn:
            self.nni = scipy.interpolate.NearestNDInterpolator(xy, v)

    def __call__(self, x1, y1):
        ret = self.lni(np.atleast_1d(x1), np.atleast_1d(y1))
        if self.outside_nn:
            sub = np.where(np.isnan(ret))
            if len(sub[0]) > 0:
                ret[sub] = self.nni(np.atleast_1d(x1), np.atleast_1d(y1))[sub]
        return(ret)
