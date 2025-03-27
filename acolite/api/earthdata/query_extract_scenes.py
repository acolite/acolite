## def query_extract_scenes
## post query and extract scenes and urls
## written by Quinten Vanhellemont, RBINS
## 2023-10-26
## modifications: 2024-04-15 (QV) added JSON options
##                2024-04-28 (QV) added as acolite function
##                2025-03-27 (QV) scene naming based on url if producer_granule_id is not present

def query_extract_scenes(query_url, verbosity = 0, link_base = 'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/'):
    import requests, time, json, os
    from bs4 import BeautifulSoup

    urls = []
    files = []

    ## post query
    if verbosity > 0: print('Posting query {}'.format(query_url))
    response = requests.get(query_url)

    type = None
    if response.ok:
        if type is None:
            try:
                dct = json.loads(response.text)
                type = 'json'
            except:
                pass

        if type is None:
            try:
                soup = BeautifulSoup(response.text, features="html.parser")
                type = 'soup'
            except:
                pass

    if type is None:
        print('Response not identified: ')
        print(response.text)
    elif type == 'json':
        if verbosity > 0:  print('Found {} entries'.format(len(dct['feed']['entry'])))
        for e in dct['feed']['entry']:
            scene = None
            if 'producer_granule_id' in e: scene = e['producer_granule_id']
            for link in e['links']:
                if 'href' not in link: continue
                if (link['href'][0:4] == 'http') &\
                   ((link['href'][-3:] in ['.nc', '.h5']) | (link['href'].lower().endswith('.zip'))):
                    urls.append('{}'.format(link['href']))
                    ## append basename if scene is not set
                    if scene is not None:
                        files.append('{}'.format(scene))
                    else:
                        files.append(os.path.basename(link['href']))
    elif type == 'soup':

        ## find reference links - really only next is important
        self = soup.find('link', href=True, rel='self')
        if self is not None: self = self['href']
        next = soup.find('link', href=True, rel='next')
        if next is not None: next = next['href']
        last = soup.find('link', href=True, rel='last')
        if last is not None: last = last['href']


        ## run through entries
        entries = soup.find_all('entry')
        if verbosity > 0:  print('Found {} entries'.format(len(entries)))
        for a in entries:
            scene = a.find('echo:producergranuleid')
            if scene is None: continue
            print(scene)

            #url = a.find('link', rel="enclosure", type="application/x-netcdf")
            url = None
            for l in a.find_all('link', rel="enclosure", type="application/x-netcdf"):
                if verbosity > 2: print(l)
                if link_base is not None:
                    if link_base in l['href']: url = '{}'.format(l['href'])
                else:
                    if (l['href'][0:4] == 'http') &\
                       ((l['href'][-3:] in ['.nc', '.h5']) | (link['href'].lower().endswith('.zip'))):
                        url = '{}'.format(l['href'])

            ## try without x-netcdf requirement
            if url is None:
                for l in a.find_all('link', rel="enclosure"):
                    if verbosity > 2: print(l)
                    if link_base is not None:
                        if link_base in l['href']: url = '{}'.format(l['href'])
                    else:
                        if (l['href'][0:4] == 'http') &\
                           ((l['href'][-3:] in ['.nc', '.h5']) | (link['href'].lower().endswith('.zip'))):
                            url = '{}'.format(l['href'])
            if url is None: continue

            urls.append('{}'.format(l['href']))
            files.append('{}'.format(scene.text))
        if verbosity > 0: print('Extracted {} urls'.format(len(urls)))

        if next is not None:
            time.sleep(0.5)
            urls_, files_ = query_extract_scenes(next)
            if len(urls_) > 0:
                urls += urls_
                files += files_

    return(urls, files)
