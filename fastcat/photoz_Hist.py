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
        self._normalize()

    def _normalize(self):
        """
        Normalize the pdfs in self.dataset by the sum.
        The true normalization should also include a multiplication
        by np.diff(self.dz)[0], assuming constant spacing as is the case.
        But this added multiplication is not needed if it is discarded in
        cPofZ and PofZ as well
        """
        integrals = self.dataset[:, 3:].sum(axis=1)
        print integrals.shape
        #integrals *= np.diff(self.dz)[0]
        self.dataset[:, 3:] = (np.where(integrals!=0., self.dataset[:, 3:].T/integrals, 0.)).T 

        
    def readfile(self, filepath):
        """
        Reads a specifically formatted file
        This could be overloaded by inheriting classes
        """
        #read comments to get the proper formatting
        #we could also do sed -n 1,4p ../test/pzdist.txt
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
        dataset.attrs['type'] = self.typestr
        dataset.attrs['file'] = self.file

    @staticmethod
    def readH5 (dataset):
        if dataset.attrs['type']==PhotoZHist.typestr:
            return PhotoZHist(dataset.attrs['file'])
        else:
            return None
        
    def applyPhotoZ (self,arr):
        #true z
        zarr = arr['z']
        #sample mag and type randomly for now
        magbins = np.random.randint(low=len(self.mag), size=len(arr))
        typebins = np.random.randint(low=len(self.type), size=len(arr))
        #add the undefined mag and types to the array,
        #and the zbins where the true z fall
        zbins = np.searchsorted(self.z, zarr)

        indices = self.tup2id(zbins, typebins, magbins)
        arr = recfunctions.append_fields(arr,
                                         ('iz','itype', 'imag'),
                                         (zbins, typebins, magbins)
                                          )
        
        return arr

    def tup2id(self, zbins, typebins, magbins):
        nmag = len(self.mag)
        ntype = len(self.type)
        return magbins + nmag*(ntype*zbins+typebins)

    def getpdf(self, arr):
        indices = self.tup2id(arr['iz'], arr['itype'], arr['imag'])
        return self.dataset[indices, 3:]
        
    def drawPhotoZ(self, arr):
        nmag = len(self.mag)
        ntype = len(self.type)
        
        #find the true z bins corresponding to the input z
        #first attempt : no mapreduce on the file bins
        #so the for loop is suboptimal
        photoz=[]
        indices = self.tup2id(arr['iz'], arr['itype'], arr['imag'])
        for idx, ztrue in zip(indices, arr['z']):
            #remove the 3 initial indices, which are iz itype imag for the corresponding line
            # how better would a direct line extraction be, compared to loading the file as
            # self.dataset and keeping it in memory : sed -n idx+5p ../test/pzdist.txt  ?
            photoz_pdf = self.dataset[idx][3:] 
            if np.all(photoz_pdf==0.):
                return -1
            cumsum = np.cumsum( photoz_pdf )
            u = np.random.random()*cumsum[-1]
            izphot = np.searchsorted(cumsum, u)
            photoz.append( self.dz[izphot] + ztrue )

        return photoz

    def PofZ(self,arr,z,dz):
        #first implementation. Need to remove the for loop eventually
        #and to cache the pdf, possibly.
        #the integration scheme is very rough: just sum of bins.
        results=[]
        for ztrue, pdf in zip(arr['z'], self.getpdf(arr)):
            xarr = ztrue + self.dz
            yarr = pdf
            masked_pdf = np.where(np.abs(xarr-z)<dz/2, yarr, 0)
            results.append(np.sum(masked_pdf))
        return np.asarray(results) #* np.diff(self.dz)[0]

    def cPofZ(self,arr,zx):
        #first implementation. Need to remove the for loop eventually
        #and to cache the pdf, possibly.
        #the integration scheme is very rough: just sum of bins
        results=[]
        pdfs = self.getpdf(arr)
        for i, ztrue in enumerate(arr['z']):
            xarr = ztrue + self.dz
            yarr = pdfs[i]
            masked_pdf = np.where(xarr<zx, yarr, 0)
            results.append(np.sum(masked_pdf))
        return np.asarray(results) #* np.diff(self.dz)[0]

