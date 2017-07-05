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
    
    def __init__(self,sigma=None,options=None,cat=None, zcent=None):
        """On input, if sigma is a number, it is interpreted as sigma=sigma_input*(1+z).
         in stuff like applyPhotoZ.  
         If cat and zcent are specified, then then neccessary fields are generated in catalog
         In that case at least arr must be specified."""
        if (cat is not None):
            if zcent is not None:
                if 'z' not in cat.data.dtype.names:
                    cat.data=recfunctions.append_fields(cat.data,'z',zcent, usemask=False)
                else:
                    cat.data['z']=zcent
            if type(sigma)==float:
                sigmar=sigma*(cat.data['z']+1)
            else:
                sigmar=sigma
                sigma=(sigma/(1+cat.data['z'])).mean()
                sigma=int(sigma*1000)/1000. ## make a nice number to output
            if 'sigma_pz' not in cat.data.dtype.names:
                cat.data=recfunctions.append_fields(cat.data,'sigma_pz',sigmar, usemask=False)
            else:
                cat.data['sigma_pz']=sigmar
        #now sigma should be a number
        if type(sigma)==type(np.ndarray):
            raise NotImplemented

        if options is not None:
            sigma=options.pz_sigma
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


    def NofZ(self,arr,zmin,zmax,dz):
        N=(zmax-zmin)/self.sigma*5
        H,zh,sh=np.histogram2d(arr["z"],arr["sigma_pz"],N)
        ## bin centrers
        zh=0.5*(zh[:-1]+zh[1:])
        sh=0.5*(sh[:-1]+sh[1:])
        zarr=np.arange(zmin,zmax,dz)
        Nz=np.zeros(len(zarr))
        for i,zc in enumerate(zh):
            for j,sig in enumerate(sh):
                if H[i,j]>0:
                    norm=1./np.sqrt(2*np.pi)/sig
                    Nz+=H[i,j]*np.exp(-(zc-zarr)**2/(2*sig**2))*norm*dz
        return zarr,Nz
        
            

    def cPofZ(self,arr,zx):
        """ Returns cumulative probability of PZ be at z<zx"""
        sigs=(zx-arr["z"])/arr["sigma_pz"]
        return 0.5*(1 + erf(sigs/np.sqrt(2)))
    
    def NameString(self):
        return "GaussPZ_"+str(self.sigma)
    
