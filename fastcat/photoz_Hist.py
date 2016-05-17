#
# Histogrammed PhotoZs : photo-z distribution provided by a
# 1D or multiD distribution, e.g. coming from simulations
#

import os, subprocess
import numpy as np
from numpy.lib import recfunctions
from photoz_Base import PhotoZBase
from scipy import interpolate, integrate

class PhotoZHist(PhotoZBase):
    """
    Histogrammed PhotoZs
    """
    typestr='hist'
    
    def __init__(self, filepath):
        if not os.path.exists(filepath):
            print 'No file %s present'%filepath
            raise
        self.file = filepath
        self.dataset = self.readfile(filepath)

    def readfile(self, filepath):
        """
        Reads a specifically formatted file
        This could be overloaded by inheriting classes
        """
        #read comments to get the proper formatting
        res = subprocess.check_output(
            ['sed', '-n', '/%s/p'%'#', filepath]
            ).split('\n')[:-1]
        for dim in res:
            n, m, M, s = dim.split()
            setattr(self, n[1:],
                    np.arange(float(m),float(M)+float(s),float(s)))

        #this takes a bit of time: improve?
        dataset = np.loadtxt(self.file)
        return dataset
    
    def writeH5(self, dataset):
        raise NotImplementedError

    @staticmethod
    def readH5 (dataset):
        raise NotImplementedError

    def applyPhotoZ (self,arr):
        #true z
        zarr = arr['z']
        #sample mag and type randomly for now
        magbins = np.random.randint(low=len(self.mag), size=len(arr))
        typebins = np.random.randint(low=len(self.type), size=len(arr))
        #add the undefined mag and types to the array,
        #and the zbins where the true z fall
        zbins = np.searchsorted(self.z, zarr)
        arr = recfunctions.append_fields(arr,
                                         ('iz','itype', 'imag'),
                                         (zbins, typebins, magbins)
                                          )
        photoz = self.drawPhotoZ(arr)
        #arr['z'] = np.asarray(photoz)
        arr = recfunctions.append_fields(arr,'zphot',photoz)
        return arr

    def drawPhotoZ(self, arr):
        nmag = len(self.mag)
        ntype = len(self.type)
        
        #find the true z bins corresponding to the input z
        #first attempt : no mapreduce on the file bins
        #so the for loop is suboptimal
        photoz=[]
        for rec in arr:
            #line index for this entry in the input file
            idx = rec['imag'] + \
              nmag*(ntype*rec['iz']+rec['itype']) + rec['imag']
            
            photoz_pdf = self.dataset[idx][3:] #remove the 3 initial indices, which are iz itype imag for the corresponding line
            if np.all(photoz_pdf==0.):
                return -1
            cumsum = np.cumsum(photoz_pdf)
            u=np.random.random()*cumsum[-1]
            izphot=np.searchsorted(cumsum, u)
            photoz.append(self.dz[izphot]+rec['z'])

        return photoz

    def PofZ(self,arr,z,dz):
        raise NotImplementedError

    def cPofZ(self,arr,zx):
        raise NotImplementedError

