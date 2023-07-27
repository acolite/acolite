## def bundle_test
## finds data and dim file for DIMAP bundle
## written by Quinten Vanhellemont, RBINS
## 2023-02-14

def bundle_test(bundle):
    import os

    dn = os.path.dirname(bundle)
    bn = os.path.basename(bundle)
    bn, ex = os.path.splitext(bn)

    if ex == '.dim':
        dimfile = bundle
        datfile = '{}/{}.data'.format(dn, bn)
    elif ex == '.data':
        dimfile = '{}/{}.dim'.format(dn, bn)
        datfile = bundle
    else:
        return()

    return(dimfile, datfile)
