## def order
## finds class to order Planet data with clip and ship and poll for status
## written by Quinten Vanhellemont, RBINS
## 2020-10-13
## modifications: 2026-09-22 (QV) integrated in acolite.api.planet

import numpy as np
import requests, json, os, time
from requests.auth import HTTPBasicAuth

class order():

    ## set up class
    def __init__(self, feature, odir, region, geojson_geometry, override = False, asset='analytic'):
        self.feature = feature
        self.region = region
        self.odir = odir
        self.geojson_geometry = geojson_geometry
        self.override = override

        self.product =  [{'item_ids': [self.feature['id']],
                          'item_type': self.feature['properties']['item_type'],
                          'product_bundle': asset}]

        print(self.product)
        ## make output file name
        self.oname = '{}.zip'.format('_'.join([self.product[0]['item_type'],
                                          self.product[0]['product_bundle'],
                                          self.product[0]['item_ids'][0],
                                          self.region]))
        self.ofile = '{}/{}'.format(self.odir, self.oname)

        # create an order request with clip and zip
        self.request_clipzip = {
              "name": "clip and zip {}".format(self.region),
              "products": self.product,
              "tools": [{"clip": {"aoi": self.geojson_geometry}}],
              "delivery": {"single_archive": True, "archive_type": "zip"}
            }

        # API Key stored as an env variable
        self.PL_API_KEY = os.getenv('PL_API_KEY')

        # set up requests to work with api
        self.auth = HTTPBasicAuth(self.PL_API_KEY, '')
        self.headers = {'content-type': 'application/json'}

        self.orders_url = 'https://api.planet.com/compute/ops/orders/v2'

    # functions for placing order and polling (adapted from Planet Notebook)
    def place_order(self):
        time.sleep(0.5+np.random.random_sample()*10)
        response = requests.post(self.orders_url, data=json.dumps(self.request_clipzip), auth=self.auth, headers=self.headers)
        if not response.ok:
            print(response.text)
            raise Exception(response.content)
        self.order_id = response.json()['id']
        self.order_url = self.orders_url + '/' + self.order_id

    ##
    def poll(self):
        time.sleep(0.5+np.random.random_sample()*2)
        r = requests.get(self.order_url, auth=self.auth)
        try:
            self.response = r.json()
        except:
            print(r)
            time.sleep(0.6)
            #self.poll()
            self.response = {'state':'timeout'}
            time.sleep(0.5+np.random.random_sample()*5)

    def poll_for_success(self, num_loops=60, sleep=30):
        time.sleep(0.5+np.random.random_sample()*2)
        count = 0
        while(count < num_loops):
            count += 1
            self.poll()
            state = self.response['state']
            print('{} - {}'.format(count,state), end='\r')
            success_states = ['success', 'partial']
            if state == 'failed':
                print(self.response)
                time.sleep(0.5+np.random.random_sample()*5)
                #raise Exception(response)
            elif state in success_states:
                break
            time.sleep(sleep)


    def process(self):
        if os.path.exists(self.ofile):
            print('We have {}'.format(self.oname))
        else:
            print('Getting {}'.format(self.oname))
            time.sleep(0.5)
            try:
                self.place_order()
            except:
                print('Could not place order')
                return()

            time.sleep(0.5)
            self.poll_for_success()


            if True:
                time.sleep(1)
                r = requests.get(self.order_url, auth=self.auth)
                response = r.json()
                state = response['state']
                time.sleep(1)
                if state != 'success':
                    print('Polling time ran out, running once more')
                    self.poll_for_success()

            ## get results from order url
            time.sleep(1)
            r = requests.get(self.order_url, auth=self.auth)
            #print(r)

            response = r.json()
            if 'results' in response['_links']:
                results = response['_links']['results']

                ## download results
                for r in results:
                    url = r['location']
                    name = r['name']

                    if '/manifest.json' in name:
                        continue

                    ## download the file
                    if self.override or not os.path.exists(self.ofile):
                        print('Downloading {} to {}'.format(name, self.ofile))
                        time.sleep(1)
                        r = requests.get(url, allow_redirects=True)
                        ## check output directory exists
                        if not os.path.exists(os.path.dirname(self.ofile)):
                            os.makedirs(os.path.dirname(self.ofile))
                        open(self.ofile, 'wb').write(r.content)
                        print('Downloaded {}'.format(self.ofile))
                    else:
                        print('{} exists'.format(self.ofile))
