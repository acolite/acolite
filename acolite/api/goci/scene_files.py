## def scene_files
## finds GOCI2 files for a given scene
## written by Quinten Vanhellemont, RBINS
## 2025-07-16
## modifications:

def scene_files(isodate, scene, dataset = 'L1B_COMP', base_url = 'https://nosc.go.kr/opendap/hyrax/GOCI-II/'):
    import requests
    from bs4 import BeautifulSoup

    ymd = isodate.split('-')

    ## set up url and get response
    url = '{}/{}/{}/{}/{}/{}/{}'.format(base_url, ymd[0], ymd[1], ymd[2], dataset, scene, 'contents.html')
    print('Accessing: {}'.format(url))
    session = requests.Session()
    response = session.get(url)
    session.close()

    ## parse response
    soup = BeautifulSoup(response.text, "html.parser")

    ## find the first table (this table contains nested tables, but we will skip those for now)
    table = soup.find('table')

    ## get table header
    header = []
    hd = table.find_all('th')
    for h in hd:
        ch = [h.text.strip() for e in h]
        header.append(ch[0])

    ## get table data
    scene_urls, scene_sizes = [], []
    rows = table.find_all('tr')
    for row in rows:
        cols_ = row.find_all('td')
        if len(cols_) == 0: continue

        ## links_ = [e.find('a', itemprop="contentUrl") for e in cols_] ## this skips the jpg files
        links_ = [e.find('a') for e in cols_] ## this includes the jpg file
        links = [l for l in links_ if l]
        if len(links) == 0: continue

        link_text = [l.text for l in links]
        link_href = [l.get('href') for l in links]

        for li, text in enumerate(link_text):
            if (text == 'file') | (text.endswith('.jpg')):
                scene_link = link_href[li]
                scene_url = '{}/{}/{}/{}/{}/{}/{}'.format(base_url, ymd[0], ymd[1], ymd[2], dataset, scene, scene_link)
                if scene_url not in scene_urls: ## with the nested table we may get doubles
                    scene_urls.append(scene_url)
                    scene_sizes.append(int(cols_[2].text.strip()))

    return(scene_urls, scene_sizes)
