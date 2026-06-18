## def rgb.basic
## function to plot 3 band scaled array
##
## written by Quinten Vanhellemont, RBINS
## 2026-06-11
## modifications: 

def basic(im, lon = None, lat = None, ofile = None, title = None, lon_lim = None, lat_lim = None,
                  axis_position = [0, 0.0, 1, 1], figsize = None, dpi = 300, facecolor = 'white'):
    import matplotlib.pyplot as plt
    import os

    ## set up plot
    fig, ax = plt.subplots(figsize = figsize)

    ## switch off first axis
    plt.axis('off')

    ## add new axes
    cax = ax.inset_axes(axis_position)

    ## make colour mesh
    if (lat is not None) & (lon is not None):
        mesh = cax.pcolormesh(lon, lat, im)

        if lon_lim is not None:
            cax.set_xlim(lon_lim)
        if lat_lim is not None:
            cax.set_ylim(lat_lim)
        cax.set_xlabel('Longitude (°E)')
        cax.set_ylabel('Latitude (°N)')
    else:
        cax.imshow(im)
        cax.axis('off')

    ## position for colourbar
    #cax = ax.inset_axes([0.93, 0.05, 0.05, 0.9])#, transform=ax.transData)
    #plt.colorbar(mesh, label = label, cax=cax, orientation = 'vertical')

    if title is not None: plt.title(title)
    if ofile is not None:
        if not os.path.exists(os.path.dirname(ofile)):
            os.makedirs(os.path.dirname(ofile))
        plt.savefig(ofile, dpi = dpi, facecolor = facecolor) #bbox_inches = 'tight'
