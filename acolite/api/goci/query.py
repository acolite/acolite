## def query
## finds and downloads GOCI2 scenes for a given date
## written by Quinten Vanhellemont, RBINS
## 2025-07-16
## modifications:

def query(date = None, scene = None, file_types = ['LA.nc'], output = None,
          download = True, override = False,
          dataset = 'L1B_COMP', base_url = 'https://nosc.go.kr/opendap/hyrax/GOCI-II/'):
    import os
    import dateutil.parser
    import acolite as ac

    ## set date if scene is given
    if scene is not None:
        scene_query = scene[0:29].replace('GOCI2', 'GC2')
        dt = scene[14:22]
        date = '-'.join((dt[0:4], dt[4:6], dt[6:8]))

    if (date is None) & (scene is None):
        print('No date or scene specified.')
        return

    ## make a list of file types to download
    if type(file_types) is not list: file_types = [file_types]

    ## download to current directory
    if output is None: output = os.getcwd()

    ## make sure we have iso date
    dt = dateutil.parser.parse(date)
    isodate = dt.isoformat()[0:10]

    ## find scenes for this date
    scene_list = ac.api.goci.date_scenes(isodate, dataset = dataset, base_url = base_url)

    ## get files matching to scene
    if scene is not None:
        if scene_query not in scene_list:
            print('Could not find scene {} for date {}'.format(scene, isodate))
        else:
            scene_list = [scene_query]

    ## get urls for each scene
    scene_base_urls = {}
    for scene_base in scene_list:
        print('Finding files for {}'.format(scene_base))
        scene_urls, scene_sizes = ac.api.goci.scene_files(date, scene_base)
        scene_base_urls[scene_base] = [z for z in zip(scene_urls, scene_sizes)]


    ## list local files
    local_files = []
    remote_files = []

    ## get files
    for scene_base in scene_base_urls:
        for dl_url, dl_size in scene_base_urls[scene_base]:
            for ft in file_types:
                if dl_url.endswith(ft):
                    bn = os.path.basename(dl_url)
                    local_file = '{}/{}'.format(output, bn)

                    if (download):
                        if os.path.exists(local_file):
                            print('Local file exists: {}'.format(local_file))
                            if (override):
                                os.remove(local_file)
                                print('Deleted local file: {}'.format(local_file))
                            else:
                                local_files.append(local_file)
                                continue

                        print('Downloading {}'.format(dl_url))
                        ret = ac.shared.download_file(dl_url, local_file)
                        if os.stat(local_file).st_size == dl_size:
                            print('Size of local file matches remote size.')
                            print('Download succesful, local file: {}'.format(local_file))
                            local_files.append(local_file)
                        else:
                            print('Size of local file does not match remote size.')
                            print('Download failed, deleting local file: {}'.format(local_file))
                            os.remove(local_file)
                    else:
                        remote_files.append(dl_url)

    if download:
        return(local_files)
    else:
        return(remote_files)
