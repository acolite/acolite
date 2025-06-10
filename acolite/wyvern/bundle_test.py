## def bundle_test
## finds json and tiff files for Wyvern data
## written by Quinten Vanhellemont, RBINS
## 2025-03-04
## modifications: 2025-06-10 (QV) updated for extracted zip files

def bundle_test(bundle):
    import os
    import acolite as ac

    scene_id = None
    dn = None
    if os.path.isdir(bundle): ## extracted zip dataset
        for f in os.listdir(bundle):
            cur_file = '{}/{}'.format(bundle, f)
            if os.path.isdir(cur_file) & (f.startswith('wyvern')): ## subfolder
                if scene_id is not None:
                    print('Multiple scenes: {}, {}'.format(scene_id, cur_file))
                scene_id = os.path.basename(f)
                for g in os.listdir(cur_file):
                    if g == '{}.json'.format(scene_id): jf = '{}/{}'.format(cur_file, g)
                    if g == '{}.tiff'.format(scene_id): file = '{}/{}'.format(cur_file, g)
            elif os.path.isfile(cur_file) & (f.startswith('wyvern')): ## no subfolder
                if scene_id is None: scene_id = f[0:46]
                if f == '{}.json'.format(scene_id): jf = '{}'.format(cur_file)
                if f == '{}.tiff'.format(scene_id): file = '{}'.format(cur_file)

    else: ## path of tiff or json is provided
        dn = os.path.dirname(bundle)
        bn = os.path.basename(bundle)
        scene_id = bn[0:46]
        jf = '{}/{}.json'.format(dn, scene_id)
        file = '{}/{}.tiff'.format(dn, scene_id)

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
