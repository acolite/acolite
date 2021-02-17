## def grid_extend
## extends sentinel2 5000m geometry datasets for interpolation
## based on method found in SNAP S2Resampling: s2tbx/s2tbx-s2msi-resampler/src/main/java/org/esa/s2tbx/s2msi/resampler/S2ResamplerUtils.java
##
## written by Quinten Vanhellemont, RBINS
## 2021-02-17
## modifications:

def grid_extend(data, ex=2, ey=2, iterations=1, crop=True):
    import numpy as np

    in_shape = data.shape
    out_shape = data.shape[0]+ey, data.shape[1]+ex

    for it in range(iterations):
        if it == 0:
            tmp_src = np.zeros(out_shape)+np.nan
            tmp_src[1:-1,1:-1]=data*1
            tmp_tar = np.zeros(out_shape)+tmp_src
        else:
            tmp_src = tmp_tar*1

        ## do y direction
        for i in range(tmp_src.shape[0]):
            for j in range(tmp_src.shape[1]):
                if np.isfinite(tmp_tar[i,j]): continue

                if (i<2):
                    s = tmp_src[i:i+3, j]
                elif (i>=2) & (i<tmp_src.shape[0]-2):
                    s = tmp_src[i-2:i+3, j]
                else:
                    s = tmp_src[i-2:i+1, j]

                l = np.where(np.isfinite(s))[0]
                if len(l) == 2:
                    if np.isnan(s[0]):
                        v = s[-2]-(s[-1]-s[-2])
                    else:
                        v =  s[1]+(s[1]-s[0])
                    tmp_tar[i,j] = v

        ## do x direction
        for i in range(tmp_src.shape[0]):
            for j in range(tmp_src.shape[1]):
                if np.isfinite(tmp_tar[i,j]): continue

                if (j<2):
                    s = tmp_src[i, j:j+3]
                elif (j>=2) & (j<tmp_src.shape[1]-2):
                    s = tmp_src[i, j-2:j+3]
                else:
                    s = tmp_src[i, j-2:j+1]

                l = np.where(np.isfinite(s))[0]
                if len(l) == 2:
                    if np.isnan(s[0]):
                        v = s[-2]-(s[-1]-s[-2])
                    else:
                        v =  s[1]+(s[1]-s[0])
                    tmp_tar[i,j] = v

    if crop: tmp_tar = tmp_tar[1:-1,1:-1]
    return(tmp_tar)
