## def calibration
## read SDBSAT1 KX10 calibration data
##
## written by Quinten Vanhellemont, RBINS
## 2023-01-03
## modifications:

def calibration(calfile):
    from xml.dom import minidom

    try:
        #xmldoc = minidom.parse(calfile)
        with open(calfile, 'r', encoding="GB2312") as f:
            xml_str = ''
            for il, line in enumerate(f.readlines()):
                xml_str+=line
                line = line.strip()
        xmldoc = minidom.parseString(xml_str)
    except:
        print('Error opening metadata file.')
        sys.exit()

    caldata = {}
    cdom = xmldoc.getElementsByTagName('MII')
    if len(cdom) > 0:
        tdom = cdom[0].getElementsByTagName('bandsIDList')
        if tdom[0].firstChild is not None:
            bandlist = tdom[0].firstChild.nodeValue
            bandlist = [b.strip() for b in bandlist.split(',')]
        for b in bandlist:
            cur = {}
            gdom = cdom[0].getElementsByTagName('RADIANCE_GAIN_BAND_{}'.format(b))
            if gdom[0].firstChild is not None:
                cur['gain'] = float(gdom[0].firstChild.nodeValue)
            gdom = cdom[0].getElementsByTagName('RADIANCE_BIAS_BAND_{}'.format(b))
            if gdom[0].firstChild is not None:
                cur['bias'] = float(gdom[0].firstChild.nodeValue)
            if len(cur) > 0: caldata[b] = cur
    return(caldata)
