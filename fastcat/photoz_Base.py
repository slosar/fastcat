#
# Basic (None) photo-z object
#
import numpy as np

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
        """ This function takes the catalog array containing
            true redshifts in "z" field and generates PDFs for 
            each galaxy in the array.
            Note that by "generating" PDFs we mean generating sufficient
            information that the PDF is completelly specified. For example,
            in case of Gaussian photo-zs, this means generating relevant means
            and variances. Such per-galaxy information should be stored back in the arr,
            which ensures it is automatically saved when dumped to HDF. Any other general
            meta-data (e.g. variance, if it is fixed for all zs, should be stored saved
            and restoed in write/readHDF
            
            For base class, there is nothing to do.
        """
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
    
    def cPofZ(self,arr,zx):
        """ Returns cumulative probability of PZ be at z<zx"""
        N=len(arr)
        cP=zeros(N)
        cP[np.where(arr["z"]<zx)]=1.0
        return cP

    def iPofZ(self,arr,zmin,zmax):
        """ Returns integrated probability int_zmin^zmax p(z) """
        return self.cPofZ(arr,zmax)-self.cPofZ(arr,zmin)
    
    def NofZ(self,arr,zmin,zmax,dz):
        """ Returns zarr, N(z) by summing probs in array arange(zmin,zmax,dz)"""
        zarr=np.arange(zmin,zmax,dz)
        Nz=np.zeros(len(zarr))
        for i,z in enumerate(zarr):
            Nz[i]=(self.PofZ(arr,z,dz)).sum()
        return zarr,Nz

    def NameString(self):
        return "TrueZ"
    
