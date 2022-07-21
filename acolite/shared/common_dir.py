## def common_dir
## find a common directory name in a list of files
## written by Quinten Vanhellemont, RBINS
## 2022-07-21
## modifications:

def common_dir(files):
    import os
    #files_paths = [os.path.abspath(f) for f in files]
    #print(files_paths)
    files_dirs = [os.path.dirname(f) for f in files]
    file_dir = None
    for fd in files_dirs:
        if len(fd) == 0: continue
        if (len(fd) == 1) & (fd[0] == os.pathsep): continue

        if all([fd in f for f in files_dirs]):
            file_dir = fd
            break
    return(file_dir)
