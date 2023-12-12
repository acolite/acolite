# def wise.bundle_test
# identifies files for wise processing
# written by Raphael Mabit, ISMER-UQAR
# 2023-12-12
#

def bundle_test(bundle):
    import os, glob
    import re
    import acolite as ac

    bn = re.sub(r"-L1.*\.pix", '', os.path.basename(bundle))
    dn = os.path.dirname(bundle)

    # WISE radiometric data
    pixhdrfile = os.path.join(dn, bn + '-L1G.pix.hdr')
    pixfile = os.path.join(dn, bn + '-L1G.pix')

    # Geo-correction Look Up tables
    gluhdrfile = os.path.join(dn, bn + '-L1A.glu.hdr')
    glufile = os.path.join(dn, bn + '-L1A.glu')

    # Navigation data: altitude, heading, pitch, roll, speed
    navlogfile = os.path.join(dn, bn + '-Navcor_sum.log')

    return(pixhdrfile, pixfile, gluhdrfile, glufile, navlogfile)