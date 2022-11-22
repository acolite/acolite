## read gdal metadata from image file
##
## written by Quinten Vanhellemont, RBINS
##
## 2022-11-22
## modifications:

def read_gdal_meta(file):
    from osgeo import gdal
    gdal.UseExceptions()

    close = False
    if type(file) == str:
        ds = gdal.Open(file)
        close = True ## close if we open file here
    elif type(file) == gdal.Dataset:
        ds = file
    else:
        print('{} not recognised'.format(file))
        return()

    try:
        md = ds.GetMetadata_Dict()
    except:
        md = {}
    if close: ds = None

    return(md)
