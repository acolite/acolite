## def srtm15plus
## downloads the srtm15plus data, returns local path
##
## SRTM15+ is provided by Tozer et al. 2019 https://doi.org/10.1029/2019EA000658
## link to dataset: https://topex.ucsd.edu/WWW_html/srtm15_plus.html
##
## function written by Quinten Vanhellemont, RBINS
## 2022-01-09
## modifications:

def srtm15plus(path=None):
    import acolite as ac
    import os

    url = 'https://topex.ucsd.edu/pub/srtm15_plus/SRTM15_V2.3.nc'
    if path is not None:
        local_file = path
    else:
        local_file = ac.config['path']+'/external/{}'.format('SRTM15_V2.3.nc')
    local_dir = os.path.dirname(local_file)

    ## download/extract
    if not os.path.exists(local_file):
        print('Downloading {} (6.5GB)'.format(local_file))
        ret = ac.shared.download_file(url, local_file, verify_ssl=False)

    ## return local path
    if os.path.exists(local_file):
        return(local_file)
    else:
        return(None)
