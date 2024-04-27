## def download
## download scenes from EarthExplorer
## written by Quinten Vanhellemont, RBINS
## 2023-09-19
## modifications: 2023-09-20 (QV) removed lxml and use HTMLparser
##                2023-11-21 (QV) added ECOSTRESS download
##                2024-04-27 (QV) moved to acolite.api

def download(entity_list, dataset_list, identifier_list, output = None,
                  extract_tar = True, remove_tar = True, override = False, verbosity = 1):
    import os, requests, json, re, netrc, time
    import acolite as ac
    #from lxml.html import fromstring
    from html.parser import HTMLParser

    ## parser to get download url
    class MyHTMLParser(HTMLParser):
        entityid = None
        dataproductid = None
        def handle_starttag(self, tag, attrs):
            dattrs = dict(attrs)
            ## Landsat data have both entity and product in secondaryDownloadButton
            if (dattrs.get('title') == 'Download Product') &\
               (dattrs.get('class') == 'btn btn-secondary secondaryDownloadButton'):
                    self.entityid = dattrs['data-entityid']
                    self.dataproductid = dattrs['data-productid']

            ## let's assume for ECOSTRESS the HDF file is the first downloadButton
            if (dattrs.get('title') == 'Download Product') &\
               (dattrs.get('class') == 'btn btn-secondary downloadButton') &\
               (self.dataproductid is None):
                    self.dataproductid = dattrs['data-productid']
            if (dattrs.get('class') == 'downloadButtons') & (self.entityid is None):
                    self.entityid = dattrs['data-entityid']
    ## end parser

    ## get output directory
    if output is None:
        cwd = os.getcwd()
        if verbosity > 0: print('No output directory given, will download to current working directory: {}'.format(cwd))
        output = '{}'.format(cwd)

    if not os.path.exists(output):
        os.makedirs(output)

    if type(entity_list) is str: entity_list = [entity_list]
    if type(dataset_list) is str: dataset_list = [dataset_list]
    if type(identifier_list) is str: identifier_list = [identifier_list]

    ## track local files
    lfiles = []

    ## run through scenes in list
    for si, scene_id in enumerate(identifier_list):
        entity_id = entity_list[si]
        dataset_id = dataset_list[si]
        scene_id = identifier_list[si]

        ## set up local and local tarfile paths
        lfile = '{}/{}'.format(output, scene_id)
        if 'ECOSTRESS' in scene_id: ## ECOSTRESS data are not tarred
            tfile = '{}/{}'.format(output, scene_id)
        else:
            tfile = '{}/{}.tar'.format(output, scene_id)

        session = None
        ## download if we don't have the scene
        if (override) | ((not os.path.exists(lfile)) & (not os.path.exists(tfile))):
            if verbosity > 1: print('Trying download of {} (entity {} dataset {})'.format(scene_id, entity_id, dataset_id))

            ## try authentication
            if verbosity > 1: print('Getting EarthExplorer access token')
            access_token, auth = ac.api.earthexplorer.auth(return_auth = True)
            if verbosity > 1: print('Got access_token {}'.format(access_token))

            ## set up session with X-Auth-Token
            session = requests.Session()
            session.headers["X-Auth-Token"] = f"{access_token}"

            ## get csrf - this seems to be needed for EROS SSO
            response = session.get(ac.config['EARTHEXPLORER_ers']+'/login')
            csrf = re.findall(r'name="csrf" value="(.+?)"', response.text)[0]

            ## log in
            if verbosity > 1: print('Logging in to EarthExplorer')
            response = session.post(ac.config['EARTHEXPLORER_ers']+'/login',
                                    data= {"username": auth[0],"password": auth[1],"csrf": csrf}, allow_redirects=True)

            if verbosity > 1: print('Got SSO cookie {}'.format(session.cookies.get("EROS_SSO_production_secure")))

            ## get scene download page
            response = session.get('{}/options/{}/{}'.format(ac.config['EARTHEXPLORER_download'], dataset_id, entity_id))

            ### get data-entityid and data-productid from Download Product button
            if verbosity > 1: print('Getting data-entityid and data-productid from Download Product button')

            ## parse html
            #tree = fromstring(response.text)
            #entityid = None
            #dataproductid = None
            #for c in tree.find_class('btn btn-secondary secondaryDownloadButton'):
            #    if c.get('title') == 'Download Product':
            #        entityid = c.get('data-entityid')
            #        dataproductid = c.get('data-productid')

            ## parse html
            parser = MyHTMLParser()
            parser.feed(response.text)
            entityid, dataproductid = parser.entityid, parser.dataproductid
            parser.close()

            if (entityid is None) | (dataproductid is None):
                print('Could not identifiy entity/dataproduct identifiers for {}'.format(scene_id))
                continue
            else:
                ## now we can create url and find download url
                url = '{}/{}/{}/EE/'.format(ac.config['EARTHEXPLORER_download'], dataproductid, entityid)

                response = session.get(url, allow_redirects=False, stream=True, timeout=1200)
                if response.ok:
                    download_url = response.json()['url']
                else:
                    print('Could not retrieve download url for {}'.format(scene_id))
                    continue

                ## try url
                print('Downloading {}'.format(scene_id))
                response = session.get(download_url, allow_redirects=False)
                ## follow redirects - not needed?
                #while response.status_code in (301, 302, 303, 307):
                #    download_url = response.headers['Location']
                #    response = session.get(download_url, allow_redirects=False)

                ## download file
                if os.path.exists(tfile): os.remove(tfile)
                dl = session.get(download_url, verify=False, allow_redirects=True)
                print('Writing file to {}'.format(tfile))
                if (dl.ok):
                    with open(tfile, 'wb') as p:
                        for chunk in dl.iter_content(chunk_size=1024*1024):
                            if chunk: # filter out keep-alive new chunks
                                p.write(chunk)
                else:
                    print('An error occurred trying to download.')
                    continue
        else:
            print('Local copy of {} exists'.format(scene_id))

        ## extract tar file
        if (extract_tar) & ('.tar' in tfile): ## if tfile == lfile then it is not tarred
            if os.path.exists(tfile):
                print('Extracting {}'.format(tfile))
                ac.shared.extract_bundle(tfile, targ_bundle=lfile)
                if not os.path.exists(lfile):
                    print('Error extracting {}'.format(tfile))
                else:
                    print('Wrote {}'.format(lfile))

            ## remove downloaded tar file after extraction
            if remove_tar:
                if os.path.exists(tfile) & os.path.exists(lfile):
                    print('Deleting {}'.format(tfile))
                    os.remove(tfile)
                    print('Deleted {}'.format(tfile))

        ## list local paths
        if extract_tar:
            if os.path.exists(lfile): lfiles.append(lfile)
        else:
            if os.path.exists(tfile): lfiles.append(tfile)
        if session is not None:
            response = session.post(ac.config['EARTHEXPLORER_ers']+'/logout')

    return(lfiles)
