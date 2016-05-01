#
# Gaussian PhotoZs
#

import numpy as np
from numpy.lib import recfunctions
from photoz_Base import PhotoZBase
from scipy.special import erf

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
        arr=recfunctions.append_fields(arr,'sigma_pz',self.sigma*(1+arr['z']), usemask=False)
        arr['z']=np.random.normal(arr['z'],arr['sigma_pz'])
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
        return arr["z"]-5*arr["sigma_pz"],arr["z"]+5*arr["sigma_pz"]

    def PofZ(self,arr,z,dz):
        """ Returns probability of PZ be at z +-dz/2"""
        norm=1./np.sqrt(2*np.pi)/arr["sigma_pz"]
        return np.exp(-(arr["z"]-z)**2/(2*arr["sigma_pz"]**2))*norm*dz

    def cPofZ(self,arr,zx):
        """ Returns cumulative probability of PZ be at z<zx"""
        sigs=(zx-arr["z"])/arr["sigma_pz"]
        return 0.5*(1 + erf(sigs/np.sqrt(2)))
    
    def NameString(self):
        return "GaussPZ_"+str(self.sigma)
    
