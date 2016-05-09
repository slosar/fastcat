##
## Class describing window function
##
import numpy as np

from photoz_Base import PhotoZBase
from photoz_Gauss import PhotoZGauss
from photoz_DoubleGauss import PhotoZDoubleGauss
from photoz_HiddenVar import PhotoZHiddenVar


## this method does not need an object    
def readPhotoZH5(dataset):
    name=dataset.attrs['type']
    toret=PhotoZBase.readH5(dataset)
    if toret is not None: return toret
    toret=PhotoZGauss.readH5(dataset)
    if toret is not None: return toret
    toret=PhotoZDoubleGauss.readH5(dataset)
    if toret is not None: return toret
    toret=PhotoZHiddenVar.readH5(dataset)
    if toret is not None: return toret
    print "Unknown PZ type!"
    stop()

    
    

    

    
