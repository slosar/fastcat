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
    """This implements mapping of z_true into f(z_true) that is multivalued
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

    For a single F that is doubled-valued in z for a given value in F,
    we get a completelly symmetric pair of Gaussians. Therefore we
    generte multiple Fs which will naturally then generate multiply sized Gaussians.

    """
    typestr='hiddenvar'
    
    def __init__(self,sigma, zcat=[0.1,0.2], zstep=[0.6,0.65]):
        """ Specifiy sigma of the main Gaussian so that error is (1+z) sigma.
             At zstep we are indistiguishable from zcat.
             zcat and zstep are arrays, so that F returns a vector rathen 
             than a vaue.

        """
        self.zcat=np.array(zcat)
        self.zstep=np.array(zstep)
        self.N=len(zcat)
        assert(len(zcat)==len(zstep))
        self.sigma=sigma
        ## if we have N measurements, need sqrt(N) less noise per one
        self.sigma1=sigma*np.sqrt(self.N)
        self.step=(np.log(1.+self.zstep)-np.log(1+self.zcat))/self.sigma
        
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
            zcat=np.array(dataset.attrs['z_cat'])
            zstep=np.array(dataset.attrs['z_step'])
            return PhotoZHiddenVar(sigma, zcat, zstep)
        else:
            return None

    def Fofz(self, z):
        F=np.outer(np.log(1+z)/self.sigma1,np.ones(self.N),)
        for i,zstep in enumerate(self.zstep):
            F[np.where(z>zstep),i]-=self.step[i]
        return F

    def dFdz(self, z):
        return 1/(1+z)*self.sigma

    
    def applyPhotoZ (self,arr,addErrors=True):
        print "Applying Hidden Var PZs"
        N=len(arr)
        F=self.Fofz(arr['z'])
        if addErrors:
            F+=np.random.normal(0,1.0,N*self.N).reshape((N,self.N))
        ## rename z into F
        nms=list(arr.dtype.names)
        nms[nms.index('z')]='Fz0'
        arr.dtype.names=tuple(nms)
        arr['Fz0']=F[:,0]
        for i in range(1,self.N):
            arr=recfunctions.append_fields(arr,'Fz'+str(i),F[:,i]) 
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
        Fm=np.zeros((len(arr),self.N))
        for i in range(self.N):
            Fm[:,i]=arr['Fz'+str(i)]
        chi2=((Ft-Fm)**2).sum(axis=1) ## error on F is one
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
        # for floats
        if (type(zx)==type(0.1)):    
            Nx=int(zx/dz)
            xar=np.linspace(0,zx,Nx)
            rect=np.zeros((Ng,Nx))
            for i,z in enumerate(xar):
                rect[:,i]=self.PofZ(arr,float(z),dz)/dz
            unnormC=trapz(rect,xar,axis=1)
        else:
        # for arrays
            zxm=zx.max()
            Nx=int(zxm/dz)
            xar=np.linspace(0,zxm,Nx)
            rect=np.zeros((Ng,Nx))
            for i,z in enumerate(xar):
                rect[:,i]=self.PofZ(arr,float(z),dz)/dz
                rect[np.where(zx>z),i]=0.0
            unnormC=trapz(rect,xar,axis=1)
        return unnormC/self.cnorm
    
    def NameString(self):
        return "HiddenVarPZ_%f_%i"%(self.sigma,self.N)
    
