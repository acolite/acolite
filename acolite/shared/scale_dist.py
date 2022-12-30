## return 'nicely rounded' scale bar lengths
##
## written by Quinten Vanhellemont, RBINS
## 2018-03-30
## modifications:

def scale_dist(scalelen):
    from math import ceil

    if scalelen < 0.1:
        sf = 1000
        md = [1, 2, 5, 10]
        unit = 'm'
    elif scalelen < 1:
        sf = 1000
        md = [100, 250, 500, 1000]
        unit = 'm'
    elif scalelen < 10:
        sf=1
        md = [1, 2, 5, 10]
        unit = 'km'
    else:
        sf=1
        md = [5, 10, 25, 50, 100]
        unit = 'km'

    scalelen*=sf
    mod = [scalelen % d for d in md]
    sel = [ri for ri,r in enumerate(mod) if r == min(mod)]
    sel = sel[-1]
    scale_nice = md[sel] * int(scalelen / md[sel])

    if scale_nice == 0: md[sel] * ceil(scalelen / md[sel])

    return(scale_nice, unit)
