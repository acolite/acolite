## def geom_grid
## makes 2D float array from sentinel metadata angles grid
## written by Quinten Vanhellemont, RBINS
## 2017-04-18
## modifications: 2021-02-17 (QV) renamed from safe_tile_grid

def grid_geom(dom):
    from numpy import zeros

    ## get column and row step
    col_step = dom.getElementsByTagName('COL_STEP')[0].firstChild.nodeValue
    row_step = dom.getElementsByTagName('ROW_STEP')[0].firstChild.nodeValue

    ## get grid values
    values = []
    for val in dom.getElementsByTagName('Values_List')[0].getElementsByTagName('VALUES'):
        values.append([float(i) for i in val.firstChild.nodeValue.split(' ')])

    ## make array of values
    nx, ny = len(values), len(values[0])
    arr =  zeros((nx,ny))
    for i in range(nx):
        for j in range(ny):
            arr[i,j] = values[i][j]
    return(arr)
