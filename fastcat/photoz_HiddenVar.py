#
# HiddenVar PhotoZ.
# This implements mapping of z_true into f(z_true) that is multivalued
# so that measurement of f(z_true) gives multiple bumps that
# by construction need to satify all probability rules

import numpy as np
from numpy.lib import recfunctions
from photoz_Base import PhotoZBase
from  scipy.integrate import trapz
from scipy.special import erf

class PhotoZHiddenVar():
    """
    This implements mapping of z_true into f(z_true) that is multivalued
    so that measurement of f(z_true) gives multiple bumps that
    by construction need to satify all probability rules

    Function is such that df/dz 
    
    Df = df/dz Dz = 1
    for Dz = (1+z) sigma_z
    
    so

    df/dz = 1/((1+z) sigma_z)
    
    so

    f=ln(1+z)/sigma_z + const

    We make a step function so that f (z_step) = f(z_cat)
    
    """
    typestr='hiddenvar'
    
    def __init__(self,sigma, zcat=0.1, zstep=0.6):
        """ Specifiy sigma of the main Gaussian so that error is (1+z) sigma.
             At zstep we are indistiguishable from zcat.
        """
        self.sigma=sigma
        self.zcat=zcat
        self.zstep=zstep
        self.step=(np.log(1+zstep)-np.log(1+zcat))/self.sigma
        
    def writeH5 (self,dataset):
        dataset.attrs['type']=self.typestr
        dataset.attrs['sigma']=self.sigma
        dataset.attrs['z_cat']=self.zcat
        dataset.attrs['z_step']=self.zstep

    @staticmethod
    def readH5 (dataset):
        """ Tries to read from H5.
            If not matched, return None
        """
        ## also use old name 
        if dataset.attrs['type']==PhotoZHiddenVar.typestr:
            sigma=float(dataset.attrs['sigma'])
            zcat=float(dataset.attrs['z_cat'])
            zstep=float(dataset.attrs['z_step'])
            return PhotoZHiddenVar(sigma, zcat, zstep)
        else:
            return None

    def Fofz(self, z):
        F=np.log(1+z)/self.sigma
        if (type(z)==type(1.0)):
            if (z>self.zstep):
                F-=self.step
        else:
            F[np.where(z>self.zstep)]-=self.step
        return F

    def dFdz(self, z):
        return 1/(1+z)*self.sigma

    
    def applyPhotoZ (self,arr,addErrors=True):
        print "Applying Hidden Var PZs"
        N=len(arr)
        F=self.Fofz(arr['z'])
        if addErrors:
            F+=np.random.normal(0,1.0,N)
        ## rename z into F
        nms=list(arr.dtype.names)
        nms[nms.index('z')]='Fz'
        arr.dtype.names=tuple(nms)
        arr['Fz']=F
        return arr
    
    def getMeanRMS (self,arr):
        """ Returns mean and sqrt variance for 
            the photoz pDF, given array. 
            Note that for assymetric PZ, you are at your
            own risk.
        """
        print "You are trying get mean rms for a gaussian error with catastrophic errorbar."
        print "I refuse to do so, but feel free to give me hard time about it"
        raise NotImplementedError

    def getMinMax(self,arr):
        """ Returns range of redshifts where p(z) is considerable, i.e.
        no real probability at z<zmin or z>zmax.
        Note that in case of catastrophic outliers, one can have considerable
        amounts of zeros in this range
        """
        print "This needs some thought to do it efficiently."
        raise NotImplementedError
        

        
    def PofZ(self,arr,z,dz):
        """ Returns probability of PZ be at z +-dz/2"""
        Ft=self.Fofz(z)
        chi2=(Ft-arr['Fz'])**2 ## error on F is one
        return np.exp(-chi2/2)*dz
    
    def cPofZ(self,arr,zx):
        ## hardcoding zmax for the time being, should fix it
        zmax=1.5
        Ng=len(arr)
        dz=0.001
        if not hasattr(self,"cnorm"):
            Nx=int(zmax/dz)
            xar=np.linspace(0,zmax,Nx)
            rect=np.zeros((Ng,Nx))
            for i,z in enumerate(xar):
                rect[:,i]=self.PofZ(arr,float(z),dz)/dz
            self.cnorm=trapz(rect,xar,axis=1)
        Nx=int(zx/dz)
        xar=np.linspace(0,zx,Nx)
        rect=np.zeros((Ng,Nx))
        for i,z in enumerate(xar):
            rect[:,i]=self.PofZ(arr,float(z),dz)/dz
        unnormC=trapz(rect,xar,axis=1)
        return unnormC/self.cnorm
    
    def NameString(self):
        return "HiddenVarPZ_%f_%f_%f"%(self.sigma,self.zcat,self.zstep)
    
