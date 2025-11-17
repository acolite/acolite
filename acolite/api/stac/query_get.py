## def query
## recursive stac query listing S2 zarr data
## written by Quinten Vanhellemont, RBINS
## 2025-11-17
## modifications:

def query_get(url, features = None):
    import requests, json
    import acolite as ac

    ## get response
    print('Querying {}'.format(url))
    response = requests.get(url)
    dict = json.loads(response.text)

    ## add features
    if 'features' in dict:
        if features is None:
            features = [f for f in dict['features']]
        else:
            features += dict['features']

    ## recurse if more pages are available
    next_url = None
    if 'links' in dict:
        for l in dict['links']:
            if l['rel'] == 'next':
                #next_url = l['href']
                extra_features = ac.api.stac.query_get(l['href'])
                for f in extra_features:
                    if f not in features:
                        features.append(f)

    ## add feature
    if 'id' in dict:
        if features is None:
            features = [dict]
        else:
            features += [dict]

    return(features)
