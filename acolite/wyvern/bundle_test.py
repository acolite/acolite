## def bundle_test
## finds json and tiff files for Wyvern data
## written by Quinten Vanhellemont, RBINS
## 2025-03-04
## modifications:

def bundle_test(bundle):
    import os
    import acolite as ac

    dn = os.path.dirname(bundle)
    bn = os.path.basename(bundle)
    scene_id = bn[0:-5]

    jf = bundle[0:-5] + '.json'
    file = bundle[0:-5] + '.tiff'

    ## download json file
    if (not os.path.exists(jf)) & False:
        base_url = 'https://wyvern-prod-public-open-data-program.s3.ca-central-1.amazonaws.com/product-type'
        if 'dragonette-001' in bn:
            product_type = 'standard'
        else:
            product_type = 'extended'
        url = '{}/{}/{}/{}.json'.format(base_url, product_type, scene_id, scene_id)
        print('Getting {}'.format(url))
        ac.shared.download_file(url, jf)

    if os.path.exists(jf) & os.path.exists(file):
        return(file, jf)
