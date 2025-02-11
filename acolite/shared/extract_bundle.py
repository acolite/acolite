## def extract_bundle
## extract a given bundle if it is zip or tar
## written by Quinten Vanhellemont, RBINS
## 2022-07-21
## modifications: 2022-10-03 (QV) check if there is only one directory in the unzipped file
##                2025-02-11 (QV) added bz2 extraction

def extract_bundle(bundle, output=None, targ_bundle=None, verbosity=0):
    import os, tarfile, zipfile, bz2, shutil
    import acolite as ac

    orig_bundle = os.path.abspath(bundle)
    bn, ex = os.path.splitext(os.path.basename(orig_bundle))
    if output is None: output = os.path.dirname(orig_bundle)
    targ_bundle = None
    extracted_path = False
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
            if verbosity>1: print("Extraction error {}, {}".format(err, type(err)))
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
                    print('Extracting {} files from {} to {}'.format(len(files), orig_bundle, targ_bundle))
                n_extracted = 0
                for z in files:
                    if verbosity>2: print('Extracting {} to {}'.format(z, targ_bundle))
                    f.extract(z, targ_bundle)
                    n_extracted += 1
            if verbosity>0: print('Extraction of {} files completed to {}'.format(n_extracted, targ_bundle))
            if file_dir is not None: targ_bundle = '{}/{}'.format(output, file_dir)
            extracted = True
        except BaseException as err:
            if verbosity>1: print("Extraction error {}, {}".format(err, type(err)))
            pass
    ## end zip file

    ## bz2 file
    if not extracted:
        try:
            targ_bundle = '{}/{}'.format(output, bn)
            with bz2.BZ2File(orig_bundle) as fi, open(targ_bundle,"wb") as fo:
                shutil.copyfileobj(fi,fo)
            if verbosity>0: print('Extraction of 1 file completed to {}'.format(targ_bundle))
            extracted = True
        except BaseException as err:
            if verbosity>1: print("Extraction error {}, {}".format(err, type(err)))
            pass
    ## end bz2 file

    if targ_bundle != None:
        ## if there is only one directory in the unzipped file
        ## use it as targ_bundle, track extracted path separately
        extracted_path = os.path.abspath(targ_bundle)
        if os.path.isdir(targ_bundle):
            dfiles = [f for f in os.listdir(targ_bundle) if f not in ['.DS_Store']]
            if len(dfiles) == 1:
                fpath = os.sep.join([targ_bundle, dfiles[0]])
                if os.path.isdir(fpath): targ_bundle = fpath

    return(targ_bundle, extracted_path)
