##
## Class describing window function
##
import numpy as np
import astropy.coordinates as ac
import healpy as hp

from window_Base import WindowBase
from window_DecBcut import WindowDecBcut
from window_Healpix import WindowHealpix
from window_Humna import WindowHumna

## this method does not need an object    
def readWindowH5(dataset):
    ## note that humna map is just a wrapper and gets saved/realoaded from healpix map
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


def registerOptions(parser):
    parser.add_option("--wf_type",dest="wftype",type="string",
                help="window func type [none,radecbcut,healpix, humna]",
                  default="humna")
    WindowDecBcut.registerOptions(parser)
    WindowHealpix.registerOptions(parser)
    WindowHumna.registerOptions(parser)

    
def getWindowFunc(o):
    if o.wftype=="none":
        return WindowBase()
    elif o.wftype=="radecbcut":
        return WindowDecBcut(options=o)
    elif o.wftype=="healpix":
        return WindowHealpix(options=o)
    elif o.wftype=="humna":
        return WindowHumna(options=o)
    else:
        print "Bad WF type:",o.wftype
        stop()
    
