#
# Gaussian PhotoZs
#

import numpy as np
from numpy.lib import recfunctions
from photoz_Gauss import PhotoZGauss
from scipy.special import erf

class PhotoZTwoPop(PhotoZGauss):
    """
    Idealised Gaussian PhotoZs, but having two popoulations,
    some with considerably worse photo-zs.
    """
    typestr='twopop'

    @staticmethod
    def registerOptions (parser):
        parser.add_option("--pz_fbad", dest="fbad", default=0.75,
                          help="PZ: fracion poor PZ for twopop", type="float")
        parser.add_option("--pz_facbad", dest="facbad", default=5.0,
                          help="PZ: factor poor PZ for twopop", type="float")
    
    def __init__(self,sigma=None,fbad=None, facbat=None, options=None):
        """ 
        A fbad fraction of galaxies will have sigma that is worse
        by facbad
        """

        if options is not None:
            sigma,fbad,facbad,=options.pz_sigma, options.fbad, options.facbad
        self.fbad=fbad
        self.facbad=facbad
        self.sigma=sigma

    def writeH5 (self,dataset):
        dataset.attrs['type']=self.typestr
        dataset.attrs['sigma']=self.sigma
        dataset.attrs['fbad']=self.fbad
        dataset.attrs['facbad']=self.facbad


    @staticmethod
    def readH5 (dataset):
        """ Tries to read from H5.
            If not matched, return None
        """
        ## also use old name 
        if (dataset.attrs['type']==PhotoZGauss.typestr):
            sigma=float(dataset.attrs['sigma'])
            fbad=float(dataset.attrs['fbad'])
            facbad=float(dataset.attrs['facbad'])
            return PhotoZTwoPop(sigma,fbad,facbad)
        else:
            return None
        
    def applyPhotoZ (self,arr):
        print "Applying PZs"
        N=len(arr)
        sigarr=self.sigma*(1+arr['z'])
        bads=np.array((np.random.uniform(0,1.,N)<self.fbad),dtype=float)
        sigarr=self.sigma*(1+arr['z'])*(1+bads*self.facbad)
        arr=recfunctions.append_fields(arr,'sigma_pz',sigarr, usemask=False)
        arr['z']=np.random.normal(arr['z'],arr['sigma_pz'])
        return arr
            
    def NameString(self):
        return "TwoPopPZ_"+str(self.sigma)+"_"+str(self.fbad)+"_"+str(self.facbad)
    
