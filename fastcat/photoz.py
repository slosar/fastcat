##
## Class describing window function
##
import numpy as np

## this method does not need an object    
def readPhotoZH5(dataset):
    name=dataset.attrs['type']
    toret=PhotoZBase.readH5(dataset)
    if toret is not None: return toret
    toret=PhotoZGauss.readH5(dataset)
    if toret is not None: return toret
    print "Unknown PZ type!"
    stop()


class PhotoZBase(object):
    """
    Basic object for describing PhotoZ error distribution.
    
    Eventually we'll probably need at least two methods to 
    work out best universal interface (one via full P(z), other via
    some sort "PZ bin")

    Base object corresponds to perfectly known PZ.

    """
    def __init__(self):
        self.type="base"

    def writeH5 (self,dataset):
        dataset.attrs['type']=self.type

    @staticmethod
    def readH5 (dataset):
        """ Tries to read from H5.
            If not matched, return None
        """
        if dataset.attrs['type']=="base":
            return PhotoZBase()
        else:
            return None

    def applyPhotoZ (self,arr):
        """ nothing to do for base class"""
        pass
        return arr

class PhotoZGauss(PhotoZBase):
    """
    Idealised Gaussian PhotoZs
    """
    typestr='gauss'
    
    def __init__(self,sigma):
        self.sigma=sigma

    def writeH5 (self,dataset):
        dataset.attrs['type']=self.typestr
        dataset.attrs['sigma']=self.sigma

    @staticmethod
    def readH5 (dataset):
        """ Tries to read from H5.
            If not matched, return None
        """
        ## also use old name 
        if (dataset.attrs['type'] in [PhotoZGauss.typestr,"fixed_sigma"]):
            sigma=float(dataset.attrs['sigma'])
            return PhotoZGauss(sigma)
        else:
            return None
        
    def applyPhotoZ (self,arr):
        print "Applying PZs"
        N=len(arr)
        arr['z']=arr['z']+(1.+arr['z'])*np.random.normal(0.,self.sigma)
        return arr
        
