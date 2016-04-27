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

    def getMeanRMS (self,arr):
        """ Returns mean and sqrt variance for 
            the photoz pDF, given array
            Note that for assymetric PZ, you are at your
            own risk.
        """
        # in base class we return redshift and zero varinace.
        N=len(arr)
        return arr["z"],np.zeros(N)

    def getMinMax(self,arr):
        """ Returns range of redshifts where p(z) is considerable, i.e.
        no real probability at z<zmin or z>zmax.
        Note that in case of catastrophic outliers, one can have considerable
        amounts of zeros in this range
        """
        return arr["z"], arr["z"]

    def PofZ(self,arr,z,dz):
        """ Returns probability of PZ be at z +-dz/2"""
        N=len(arr)
        P=zeros(N)
        dzo2=dz/2
        P[np.where(abs(arr["z"]-z)<dzo2)]=1.0
        return P

    def NameString(self):
        return "TrueZ"
    
    
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
        
    def getMeanRMS (self,arr):
        """ Returns mean and sqrt variance for 
            the photoz pDF, given array. 
            Note that for assymetric PZ, you are at your
            own risk.
        """
        return arr["z"],self.sigma*(1+arr["z"])

    def getMinMax(self,arr):
        """ Returns range of redshifts where p(z) is considerable, i.e.
        no real probability at z<zmin or z>zmax.
        Note that in case of catastrophic outliers, one can have considerable
        amounts of zeros in this range
        """
        return arr["z"]-self.sigma*(1+arr["z"])*5,arr["z"]+self.sigma*(1+arr["z"])*5

    def PofZ(self,arr,z,dz):
        """ Returns probability of PZ be at z +-dz/2"""
        sig=self.sigma*(1+arr["z"])
        norm=1./np.sqrt(2*np.pi)/sig
        return np.exp(-(arr["z"]-z)**2/(2*self.sigma**2))*norm*dz

    def NameString(self):
        return "GaussPZ_"+str(self.sigma)
    
    
