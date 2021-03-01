## def settings_write
## write "settings" file for ACOLITE based on settings dict
##
## Written by Quinten Vanhellemont 2017-11-30
## Last modifications:

def write(file, settings):
    import os, datetime
    if not os.path.exists(os.path.dirname(file)): os.makedirs(os.path.dirname(file))
    comf='## {}\n'
    valf='{}={}\n'
    with open(file,'w') as f:
        f.write(comf.format('ACOLITE settings'))
        f.write(comf.format('Written at {}'.format(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))))
        for key in settings.keys():
            vals = settings[key]
            if type(vals) is list:
                vals = ','.join([str(i) for i in settings[key]])
            f.write(valf.format(key,vals))
