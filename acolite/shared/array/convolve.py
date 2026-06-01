## def convolve
## convolution of an array across a given axis
## written by Quinten Vanhellemont, RBINS
## 2026-06-01
## modifications:

def convolve(x, window = 5, axis = 0, iterations = 3):
    import scipy.signal
    import numpy as np

    ## set up kernel
    kernel_dim = [1] * len(x.shape)
    kernel_dim[axis] = window
    kernel = np.ones(kernel_dim)/window

    ## convolute using fft
    x_ = x * 1.0
    for it in range(iterations):
        x_ = scipy.signal.fftconvolve(x_, kernel, 'same')
    return(x_)
