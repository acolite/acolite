## def closest_idx
## returns closest index and selected value from list
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-08
## modifications:
##                

def closest_idx(xlist, xval):
    idx,xret = min(enumerate(xlist), key=lambda x: abs(float(x[1])-float(xval)))
    return(idx,xret)
