## def bundle tests
## finds image files and metadata files for AMAZONIA bundle
## written by Quinten Vanhellemont, RBINS
## 2022-01-14

def bundle_test(bundle):
    import glob
    files_tiff = glob.glob('{}/*.tif'.format(bundle))
    files_tiff.sort()
    files_xml = glob.glob('{}/*.xml'.format(bundle))
    files_xml.sort()

    return(files_xml, files_tiff)
