## def bundle_test
## finds if given directory is a Pléiades image bundle
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-06-29
## modifications: QV 2017-03-16 added nc metadata
##                QV 2017-03-23 new parser & changed method to support multiple scenes in same package
##                QV 2017-05-04 nc returns none lists
##                QV 2019-02-27 check pan file for tile
##                QV 2021-02-24 renamed from image_test removed nc option


def bundle_test(file, listpan=True):
    import os, sys, fnmatch
    from xml.dom import minidom

    sensor=[]
    product=set()
    img_dir=[]
    xml_file=[]
    img_file=[]
    sub_file = []

    ## bundle metadata
    for i, fname in enumerate(os.listdir(file)):
        fname_ = fname.upper()

        ## bundle metadata
        if fnmatch.fnmatch(fname_,'VOL_PHR.XML') or fnmatch.fnmatch(fname_,'SPOT_VOL.XML'):
            xmldoc = minidom.parse('{}/{}'.format(file,fname))
            components=[]
            for t in xmldoc.getElementsByTagName('Component'):
                node = t.getElementsByTagName('COMPONENT_TITLE')
                if len(node) > 0: title = node[0].firstChild.nodeValue
                node = t.getElementsByTagName('COMPONENT_PATH')
                if len(node) > 0: path = node[0].attributes['href'].value
                tmp = title.split(' ')
                if tmp[0] not in ['SENSOR','ORTHO','SYSTEM_ORTHO']: continue
                _, dtype, dataset = tmp
                components.append({'title':title,'path':path, 'type': dtype, 'dataset': dataset})

        ## image directories
        if fnmatch.fnmatch(fname_,'IMG*'):
            img = file+'/'+fname
            sub_file.append(fname)
            img_dir.append(img)
            split = fname.split('_')
            sensor.append(split[1])
            product.add(split[2])
            for j, ffile in enumerate(os.listdir(img)):
                ffile_ = ffile.upper()
                if fnmatch.fnmatch(ffile_, 'DIM_'+'*.XML'): xml_file.append(img+'/'+ffile)
                if fnmatch.fnmatch(ffile_, 'IMG_'+'*.TIF'): img_file.append(img+'/'+ffile)
                if fnmatch.fnmatch(ffile_, 'IMG_'+'*.JP2'): img_file.append(img+'/'+ffile)

    if len(product) > 1:
        if ('P' in product) and ('MS' in product): ptype = 'bundle'
    else:
        if ('PMS' in product): ptype = 'PMS'
        if ('MS' in product): ptype = 'MS'

    imagefile = []
    metafile= []
    imagedataset = []
    panimagefile = []
    panmetafile = []

    if (ptype == 'bundle') | (ptype == 'MS'):
        for i, img in enumerate(img_file):
            sel_image = ''
            sel_meta = ''
            sel_dataset = ''
            sel_panimage = ''
            sel_panmeta = ''

            if (fnmatch.fnmatch(img.upper(), '*IMG*MS*')):
                sf = img.split('/')[-2]
                md = [md for md in xml_file if (sf in md) & (fnmatch.fnmatch(md.upper(), '*IMG*MS*.XML'))]

                ## we have metadata and image
                if (len(md) == 1):
                    sel_image = img # imagefile.append(im[0])
                    sel_meta = md[0] # metafile.append(md[0])
                    for component in components:
                        if sf == component['path'].split('/')[0]:
                            sel_dataset=component['dataset']

                ## find tile to determine proper pan file (this is not perfect)
                bn = os.path.basename(img)
                bn = bn[0:bn.find('.')]
                tile = bn.split('_')[-1]

                ## find matching pan
                for component in components:
                    if (component['dataset'] == sel_dataset) & (component['type'] == 'P'):
                        psf = component['path'].split('/')[0]
                        for pi, pimg in enumerate(img_file):
                            if (psf in pimg) & (fnmatch.fnmatch(pimg.upper(), '*IMG*_P_*_{}.*'.format(tile))):
                                pmd = [md for md in xml_file if (psf in md) & (fnmatch.fnmatch(md.upper(), '*IMG*_P_*.XML'))]
                                if (len(pmd) == 1):
                                    sel_panimage=pimg
                                    sel_panmeta=pmd[0]

                imagefile.append(sel_image)
                metafile.append(sel_meta)
                imagedataset.append(sel_dataset)
                panimagefile.append(sel_panimage)
                panmetafile.append(sel_panmeta)

    if (ptype == 'PMS'):
        for i, sf in enumerate(sub_file):
            im = [im for im in img_file if (sf in im) & (fnmatch.fnmatch(im.upper(), '*IMG*PMS*'))]
            md = [md for md in xml_file if (sf in md) & (fnmatch.fnmatch(md.upper(), '*IMG*PMS*.XML'))]
            if (len(im) == 1) & (len(md) == 1):
                imagefile.append(im[0])
                metafile.append(md[0])
                panimagefile.append('')
                panmetafile.append('')

    if (len(metafile) != 0) & (len(imagefile) != 0):
        if listpan:
             return(imagefile, metafile, panimagefile, panmetafile)
        else:
             return(imagefile, metafile)
    else:
        print('Directory '+file+' probably not a Pléiades image.')
        return()
