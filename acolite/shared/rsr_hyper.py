## def rsr_hyper
## makes hyperspectral RSR based on center wavelength and band width
## written by Quinten Vanhellemont, RBINS
## 2021-06-08
## modifications: 2024-07-03 (QV) added type and factor keywords
##                2025-11-26 (QV) added names keyword

def rsr_hyper(waves, widths, step = 0.25, nm = True, type = 'gauss', factor = None, names = None):
    import acolite as ac

    nbands = len(waves)
    rsr = {}
    for b in range(nbands):
        if type == 'gauss':
            if factor is None: factor = 1.5
            wave, resp = ac.shared.gauss_response(waves[b], widths[b], step = step, factor = factor)
        elif type == 'square':
            if factor is None: factor = 0.5
            wave, resp = ac.shared.square_response(waves[b], widths[b], step = step, factor = factor)

        ## set name if names are given
        bname = b
        if names is not None:
            if len(names) == nbands: bname = names[b]

        rsr[bname] = {'wave':wave, 'response': resp}
        if nm: rsr[bname]['wave']/=1000
    return(rsr)
