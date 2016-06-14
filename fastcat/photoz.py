##
## Class describing window function
##
import numpy as np

from photoz_Base import PhotoZBase
from photoz_Gauss import PhotoZGauss
from photoz_TwoPop import PhotoZTwoPop
from photoz_HiddenVar import PhotoZHiddenVar
from photoz_Hist import PhotoZHist


## this method does not need an object    
def readPhotoZH5(dataset):
    name=dataset.attrs['type']
    toret=PhotoZBase.readH5(dataset)
    if toret is not None: return toret
    toret=PhotoZGauss.readH5(dataset)
    if toret is not None: return toret
    toret=PhotoZTwoPop.readH5(dataset)
    if toret is not None: return toret
    toret=PhotoZHiddenVar.readH5(dataset)
    if toret is not None: return toret
    toret=PhotoZHist.readH5(dataset)
    if toret is not None: return toret
    print "Unknown PZ type!"
    stop()

def registerOptions(parser):
    parser.add_option("--pz_type",dest="pztype",type="string",
        help="photo z type [none,gauss, twopop, hiddenvar, franzona]",
                  default="gauss")
    ## commonly used, so register it globally
    parser.add_option("--pz_sigma", dest="pz_sigma", default=0.01,
                      help="PZ: Guass sigma for (1+z)", type="float")
    PhotoZTwoPop.registerOptions(parser)
    PhotoZHiddenVar.registerOptions(parser)
    PhotoZHist.registerOptions(parser)
    
def getPhotoZ(o):
    if o.pztype=="none":
        pz = PhotoZBase()
    elif o.pztype=="gauss":
        pz = PhotoZGauss(options=o)
    elif o.pztype=="twopop":
        pz = PhotoZTwoPop(options=o)
    elif o.pztype=="hiddenvar":
        pz = PhotoZHiddenVar(options=o)
    elif o.pztype=="franzona":
        pz = PhotoZHist(options=o)
    else:
        print "Bad PZ type:",o.pztype
        stop()
    return pz

