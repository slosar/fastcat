##
## Class describing window function
##
import numpy as np


class WindowBase(object):
    """
    Basic object for describing window function of survey.
    It needs to implement one function: 

    float __call__(ra,dec),

    which returns relative angular probability of seeing
    an object (between 0 and 1).

    ras and decs can also be numpy arrays in which case
    it returns a numpyu array
    """
    typestr="base"
    def __init__(self):
        pass
        
    def __call__(self,ra,dec):
        if (type(ra)==type(1.0)):
            return 1.0
        else:
            return np.ones(len(ra))

    def writeH5 (self,of):
        dset=of.create_dataset("window",data=[])
        dset.attrs['type']=self.typestr

    @staticmethod
    def readH5 (dataset):
        """ Tries to read from H5.
            If not matched, return None
        """
        if dataset.attrs['type']==WindowBase.typestr:
            return WindowBase()
        else:
            return None

    def applyWindow (self,arr):
        """ Applies window function to
            a catalog's array by proportionally (randomly) cutting
            entries given the window efficiency.
        """
        ra=arr['ra']
        dec=arr['dec']
        N=len(ra)
        wf=self(ra,dec)
        mask=(wf>np.random.uniform(0,1,N))
        arr=arr[mask]
        print "After window: ",N,"->",len(arr)
        return arr

    def NameString (self):
        return "FullSky"
    
