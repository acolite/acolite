import os, sys, datetime
## object for logging stdout to log file when processing
class LogTee(object):

        def __init__(self, name):
            self.name=name
            ## make new file
            if os.path.exists(os.path.dirname(self.name)) is False:
                try:
                    os.makedirs(os.path.dirname(self.name))
                except:
                    print('Error: could not create directory: {}'.format(os.path.dirname(self.name)))
                    exit(1)
            self.file = open(self.name, 'w')
            self.file.close()
            self.mode='a'
            self.stdout = sys.stdout
            sys.stdout = self
        def __del__(self):
            try:
                sys.stdout = self.stdout
            except:
                pass
        def write(self, data):
            self.stdout.write(data)
            data = data.strip()
            if len(data) > 0:
                with open(self.name, self.mode) as self.file:
                    self.file.write('{}: {}\n'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),data))
        def flush(self):
            pass
