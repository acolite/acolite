## def date_scenes
## finds GOCI2 scenes for a given date
## written by Quinten Vanhellemont, RBINS
## 2025-07-16
## modifications:

def date_scenes(isodate, dataset = 'L1B_COMP', base_url = 'https://nosc.go.kr/opendap/hyrax/GOCI-II/'):
    import requests
    from bs4 import BeautifulSoup

    ymd = isodate.split('-')

    ## set up url and get response
    url = '{}/{}/{}/{}/{}'.format(base_url, ymd[0], ymd[1], ymd[2], dataset)
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
    scenes = []
    rows = table.find_all('tr')
    for row in rows:
        cols_ = row.find_all('td')
        if len(cols_) == 0: continue
        #cols = [e.text.strip() for e in cols_]
        ## these are from a nested subtable - skip them
        #if cols[0] in ['', '-']: continue
        #scene = cols[0].strip('/')
        #scenes.append(scene)
        ## find links
        links_ = [e.find('a') for e in cols_]
        links = [l for l in links_ if l]
        if len(links) == 0: continue
        text_ = [l.text for l in links]
        scene = text_[0].strip('/')
        scenes.append(scene)

    return(scenes)
