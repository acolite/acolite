## def bundle_test
## bundle test for DEIMOS2/GEOSAT2 data
## returns dict with paths to metadata file and iamge files
## written by Quinten Vanhellemont, RBINS
## 2024-05-16
## modifications: 2025-05-15 (QV) added Deimos1

def bundle_test(bundle):

    import os, glob
    files = glob.glob(bundle+'/*')
    files.sort()

    data = {}
    for file in files:
        bn = os.path.basename(file)

        if bn[0:11] == 'DE2_MS4_L1C':
            sensor = 'DEIMOS2_HIRAIS'
            sensor = 'GEOSAT2_HIRAIS'
            datatype = 'MS'
            level = 'L1C'

        elif bn[0:11] == 'DE2_PAN_L1C':
            sensor = 'DEIMOS2_HIRAIS'
            sensor = 'GEOSAT2_HIRAIS'
            datatype = 'PAN'
            level = 'L1C'
        elif bn[0:12] == 'DE01_SL6_22S':
            sensor = 'DEIMOS1_SLIM6'
            sensor = 'GEOSAT1_SLIM6'
            datatype = 'MS'
            level = 'L1C'
        else:
            #print(bn[0:11])
            continue

        ext = os.path.splitext(bn)[1]
        if datatype not in data: data[datatype] = {}

        if ext == '.dim': filetype = 'metadata'
        elif ext == '.tif': filetype = 'image'
        elif ext == '.png': filetype = 'quicklook'
        elif ext == '.html': filetype = 'order'
        else: filetype = ext[1:]

        if filetype in data[datatype]:
            print('Replacing', datatype, filetype, data[datatype][filetype])
            print('with', datatype, filetype, file)

        data[datatype][filetype] = file

    if len(data) == 0: return
    return(data)
