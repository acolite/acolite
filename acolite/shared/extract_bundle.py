## def extract_bundle
## extract a given bundle if it is zip or tar
## written by Quinten Vanhellemont, RBINS
## 2022-07-21
## modifications:

def extract_bundle(bundle, output=None, targ_bundle=None, verbosity=0):
    import os, tarfile, zipfile
    import acolite as ac

    bn, ex = os.path.splitext(os.path.basename(bundle))
    orig_bundle = '{}'.format(bundle)
    if output is None: output = os.path.dirname(orig_bundle)
    targ_bundle = None
    extracted = False

    ## tar file
    if not extracted:
        try:
            with tarfile.open(orig_bundle) as f:
                files = f.getnames()
                file_dir = ac.shared.common_dir(files)
                if file_dir is None:
                    targ_bundle = '{}/{}'.format(output, bn)
                else:
                    targ_bundle = '{}'.format(output)
                if verbosity>0:
                    print('Extracting {} files from {} to {}'.format(len(files), bundle, targ_bundle))
                f.extractall(targ_bundle)
            if verbosity>0: print('Extraction completed to {}'.format(targ_bundle))
            if file_dir is not None: targ_bundle = '{}/{}'.format(output, file_dir)
            extracted = True
        except BaseException as err:
            if verbosity>1: print(f"Extraction error {err=}, {type(err)=}")
            pass
    ## end tar file

    ## zip file
    if not extracted:
        try:
            with zipfile.ZipFile(orig_bundle, 'r') as f:
                files= [i.filename for i in f.infolist()]
                file_dir = ac.shared.common_dir(files)
                if file_dir is None:
                    targ_bundle = '{}/{}'.format(output, bn)
                else:
                    targ_bundle = '{}/'.format(output)
                if verbosity>0:
                    print('Extracting {} files from {} to {}'.format(len(files), bundle, targ_bundle))
                for z in files:
                    f.extract(z, targ_bundle)
            if verbosity>0: print('Extraction completed to {}'.format(targ_bundle))
            if file_dir is not None: targ_bundle = '{}/{}'.format(output, file_dir)
            extracted = True
        except BaseException as err:
            if verbosity>1: print(f"Extraction error {err=}, {type(err)=}")
            pass
    ## end zip file

    return(targ_bundle)
