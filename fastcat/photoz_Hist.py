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
        
    def join_struct_arrays(self, arrays):
        newdtype = sum((a.dtype.descr for a in arrays), [])
        newrecarray = np.empty(len(arrays[0]), dtype = newdtype)
        for a in arrays:
            for name in a.dtype.names:
                newrecarray[name] = a[name]
        return newrecarray

    def applyPhotoZ (self,arr):
        #true z
        zarr = arr['z']
        #sample mag and type randomly for now
        magbins = np.array(np.random.randint(low=len(self.mag), size=len(arr)), dtype=[('imag', np.int)])
        typebins = np.array(np.random.randint(low=len(self.type), size=len(arr)), dtype=[('itype', np.int)])
        #add the undefined mag and types to the array,
        #and the zbins where the true z fall
        zbins = np.array(np.searchsorted(self.z, zarr), dtype=[('iz', np.int)])
        newarr = self.join_struct_arrays([zbins, typebins, magbins])
        return self.join_struct_arrays([arr, newarr])

    def tup2id(self, zbins, typebins, magbins):
        nmag = len(self.mag)
        ntype = len(self.type)
        return magbins + nmag*(ntype*zbins + typebins)

    def getpdf(self, arr):
        indices = self.tup2id(arr['iz'], arr['itype'], arr['imag'])
        return self.dataset[indices, 3:]
        
    def drawPhotoZ(self, arr, nsamples):
        ngal=arr.size
        photoz_pdfs = self.getpdf(arr)
        pz_samples = np.zeros((ngal, nsamples))
        #deal immediately with flat 0 pdfs:
        mask = photoz_pdfs.sum(axis=1)!=0
        masked_pdfs = photoz_pdfs[mask]
        #pz_samples[np.logical_not(mask)] = np.zeros(shape=(nsamples,))

        cumsum = np.cumsum(masked_pdfs, axis=1)
        masked_ngal=cumsum.shape[0]
        #draw uniform random numbers
        u = np.random.random(masked_ngal*nsamples)
        u=u.reshape(nsamples, masked_ngal)*cumsum[:,-1]
        
        dz=np.zeros((masked_ngal, nsamples))
        #now loop over the masked gal to call searchsorted
        #which only accepts 1D array as input
        #this for loop will not scale well with high number of galaxies
        for i in range(masked_ngal):
            izphot = np.searchsorted(cumsum[i], u.T[i])
            dz[i,:] = self.dz[izphot]
        zphot = dz.T + arr[mask]['z']
        pz_samples[mask] = zphot.T

        return pz_samples

    def PofZ(self,arr,z,dz):
        #Need to cache the pdf?
        #the integration scheme is very rough: just sum of bins.
        results=np.zeros((arr.size))
        pdfs = self.getpdf(arr)
        mask = pdfs.sum(axis=1)!=0
        xarrs = arr[mask]['z'][:,np.newaxis] + self.dz
        masked_pdf = np.where(np.abs(xarrs-z)<dz/2., pdfs[mask], 0.)
        results[mask] = np.sum(masked_pdf, axis=1)
        return results

    def cPofZ(self,arr,zx):
        #Need to cache the pdf?
        #the integration scheme is very rough: just sum of bins
        results=np.zeros((arr.size))
        pdfs = self.getpdf(arr)
        mask = pdfs.sum(axis=1)!=0
        xarrs = arr[mask]['z'][:,np.newaxis] + self.dz
        masked_pdf = np.where(xarrs<zx, pdfs[mask], 0.)
        results[mask] = np.sum(masked_pdf, axis=1)
        return results
