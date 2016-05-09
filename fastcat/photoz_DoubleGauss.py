#
# DoubleGaussian PhotoZ.
# Still idealistic knowledge, but assymetric, bimodal, difficult photozs
#

import numpy as np
from numpy.lib import recfunctions
from photoz_Base import PhotoZBase
from scipy.special import erf

class PhotoZDoubleGauss():
    """
    Photo-zs with double Gaussian, where one of the two is fixed at a
    'catastrophic' redshift.

    This is somewhat tricky to implement right. See 

    https://docs.google.com/presentation/d/1PvDOfGqh4UT3Ulasp7KClBPzFJtF7lAeKv02I4M_igA/edit?usp=sharing
    
    """
    typestr='doublegauss'
    
    def __init__(self,sigma, Acat, zcat, sigmacat):
        """ Specifiy sigma of the main Gaussian
            and that of the catastrophic side. Acat is the probabily in the
            catastrophic Gaussian.
        """
        self.sigma=sigma
        self.Acat=Acat
        self.zcat=zcat
        self.sigmacat=sigmacat
        
    def writeH5 (self,dataset):
        dataset.attrs['type']=self.typestr
        dataset.attrs['sigma']=self.sigma
        dataset.attrs['Acat']=self.Acat
        dataset.attrs['zcat']=self.zcat
        dataset.attrs['sigmacat']=self.sigmacat

    @staticmethod
    def readH5 (dataset):
        """ Tries to read from H5.
            If not matched, return None
        """
        ## also use old name 
        if dataset.attrs['type']==PhotoZDoubleGauss.typestr:
            sigma=float(dataset.attrs['sigma'])
            Acat=float(dataset.attrs['Acat'])
            zcat=float(dataset.attrs['zcat'])
            sigmacat=float(dataset.attrs['sigmacat'])
            return PhotoZDoubleGauss(sigma,Acat,zcat,sigmacat)
        else:
            return None

    def applyPhotoZ (self,arr):
        print "Applying Double Gauss PZs"
        N=len(arr)
        Nc=np.random.binomial(N,self.Acat)
        #new z array to store stuff
        nz=np.zeros(N)
        zsigma=np.zeros(N)
        ztrue=arr['z']
        ## now sort ztrue and indices. This will be my master list of objects
        print "   creating new z array..."
        zndx=[(v[1],v[0]) for v in enumerate(ztrue)]
        print "   sorting..."
        zndx.sort()
        print "   generate cats..."
        ## first find Nc guys sampled from cat z for swapping
        catl=np.random.normal(self.zcat,self.sigmacat,Nc)
        catl.sort()
        i=0
        j=0
        upscat=[]
        ## make sure we can extend far enough
        try:
            assert(catl[0]>zndx[0][0])
            assert(catl[-1]<zndx[-1][0])
        except AssertionError:
            print "someting bad has happened: "
            print "zcat min max=",catl[0], catl[-1]
            print "scat min max=",zndx[0], zndx[-1]
            raise AssertionError
        print "   sampling cats..."
        deltaz=0
        while (i<Nc):
            while(zndx[j][0]<catl[i]):
                j+=1
            ## now take this or the one before
            j+=np.random.randint(-1,1) ## add -1 or 0
            upscat.append(zndx[j])
            deltaz=max(deltaz,abs(zndx[j][0]-catl[i]))
            i+=1
            ## step back just in case
            j-=1
            if (i%1000==0):
                print i,j,Nc,'\r',
        print "Max z error:",deltaz
        print "   actually generating cats..."
        ## now we are ready to rock and roll
        mask=np.zeros(N,dtype=bool)
        c=0
        for z,i in upscat:
            mask[i]=True
            while True:
                z2,i2=zndx[np.random.randint(0,N)]
                if not mask[i2]:
                    mask[i2]=True
                    break
            czsigma=self.sigma*(1+z2)    
            zc=np.random.normal(z2,czsigma)
            nz[i]=zc
            nz[i2]=zc
            zsigma[i]=czsigma
            zsigma[i2]=czsigma
            if (c%1000==0):
                print c,Nc,'\r',
            c+=1

        print "   actually generating noncats..."
        for z,i in zndx:
            if nz[i]==0:
                czsigma=self.sigma*(1+z2)    
                nz[i]=np.random.normal(z,czsigma)
                zsigma[i]=czsigma
                
        ## rename z in zmg
        nms=list(arr.dtype.names)
        nms[nms.index('z')]='zmg'
        arr.dtype.names=tuple(nms)
        arr['zmg']=nz
        ## add zsigma
        arr=recfunctions.append_fields(arr,'sigma_pz',zsigma, usemask=False)
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
        minz=arr['zmg']-arr['sigma_pz']*5
        dmin=self.zcat-5*self.sigmacat
        minz[np.where(minz>dmin)]=dmin
        maxz=arr['zmg']+arr['sigma_pz']*5
        dax=self.zcat+5*self.sigmacat
        maxz[np.where(maxz<dmax)]=dmax
        return dmin,dmax
        
    def PofZ(self,arr,z,dz):
        """ Returns probability of PZ be at z +-dz/2"""
        norm=(1-self.Acat)/np.sqrt(2*np.pi)/arr["sigma_pz"]
        norm2=self.Acat/np.sqrt(2*np.pi)/self.sigmacat
        return np.exp(-(arr['zmg']-z)**2/(2*arr['sigma_pz']**2))*norm+np.exp(-(self.zcat-z)**2/(2*self.sigmacat**2))*norm2

    def cPofZ(self,arr,zx):
        """ Returns cumulative probability of PZ be at z<zx"""
        sig1=(zx-arr['zmg'])/arr['sigma_pz']
        sig2=(zx-self.zcat)/self.sigmacat
        return (1.-self.Acat)*0.5*(1 + erf(sig1/np.sqrt(2)))+self.Acat*0.5*(1 + erf(sig2/np.sqrt(2)))
    
    def NameString(self):
        return "DoubleGaussPZ_%f_%f_%f_%f"%(self.sigma,self.Acat,self.zcat,self.sigmacat)
    
