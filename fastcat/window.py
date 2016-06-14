##
## Class describing window function
##
import numpy as np
import astropy.coordinates as ac
import healpy as hp

from window_Base import WindowBase
from window_DecBcut import WindowDecBcut
from window_Healpix import WindowHealpix

## this method does not need an object    
def readWindowH5(dataset):
    name=dataset.attrs['type']
    ## loop over possible types
    ## (why can't i do 'for t in [WindowBase, WindowDecBcut]:'?)
    toret=WindowBase.readH5(dataset)
    if toret is not None: return toret
    toret=WindowDecBcut.readH5(dataset)
    if toret is not None: return toret
    toret=WindowHealpix.readH5(dataset)
    if toret is not None: return toret
    print "Unknown window type!"
    stop()


