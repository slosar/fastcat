#
# Template SED photo-z object
#
import numpy as np
from numpy.lib import recfunctions
import os

import sys
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d

from photoz_Base import PhotoZBase

class PhotoZTemplateSED(PhotoZBase):
    """
    TemplateSED method
    """

    def __init__(self, infopath=None, options=None):
        self.path = infopath
        typestr='TemplateSED'
        if options is not None:
            self.sigma=options.pz_sigma
            self.f_obs=options.pz_flux

    def writeH5 (self,dataset):
        dataset.attrs['type']=self.type
        dataset.attrs['sigma']=self.sigma
        dataset.attrs['flux']=self.f_obs

    @staticmethod
    def readH5 (dataset):
        """ Tries to read from H5.
            If not matched, return None
        """
        if dataset.attrs['type']=="base":
            return PhotoZBase()
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
        print "Applying Template SED PZs"

        ztrue = arr['z']

        #select a template
        templates = ['El_B2004a.sed']+['Sbc_B2004a.sed','Scd_B2004a.sed']
        templates = templates +['Im_B2004a.sed','SB3_B2004a.sed','SB2_B2004a.sed','ssp_25Myr_z008.sed','ssp_5Myr_z008.sed']

        #read in f_mod files, interpolate, get values of f_mod_b
        self.z_grid = np.load(os.path.join(self.path,'z_grid.npy'))
        self.f_mod = np.load(os.path.join(self.path,'f_mod.npy'))
        (self.nz,self.nt,self.nb) = np.shape(self.f_mod)
        ngals = len(ztrue)

        f_mod_o = np.zeros((self.nb, ngals))
        for z in range(ngals):
            #currently templates are randomly chosen but probably should be an input with true z
            templateno = np.random.choice(range(self.nt))
            for b in range(self.nb):
                spl = InterpolatedUnivariateSpline(self.z_grid, self.f_mod[:,templateno,b])
                f_mod_o[b][z] = spl(ztrue[z])

        #select sigma_b - 10% for now
        self.sigma = 0.1*f_mod_o

        #select observed fluxes f_obs_b = f_mod_b + sigma_b*rando
        self.f_obs = np.zeros((self.nb, ngals))
        for n in range(ngals):
            for b in range(self.nb):
                self.f_obs[b][n] = f_mod_o[b][n] + self.sigma[b][n] * np.random.normal()

        #eventually want this to be meaningful output
        return self.f_obs.T,self.sigma.T

    def scalefree_flux_likelihood(self, f_obs, f_obs_var, f_mod):

        var = f_obs_var  # nz * nt * nb
        invvar = np.where(f_obs/var < 1e-6, 0.0, var**-1.0)  # not nz * nt * nb
        FOT = np.sum(f_mod * f_obs * invvar, axis=2)  # nz * nt
        FTT = np.sum(f_mod**2 * invvar, axis=2)  # nz * nt
        FOO = np.dot(invvar, f_obs**2)  # nz * nt
        chi2 = FOO - FOT**2.0 / FTT  # nz * nt
        like = np.exp(-0.5*chi2) / np.sqrt(FTT)  # nz * nt
        return like

    def getpdf(self,f_obs,sigma,z,dz):

        ff_obs = f_obs
        f_obs_err = sigma

        dzdz = dz/10.

        zrange = np.arange(max(0.,z-dz/2.),min(z+dz/2.,1.4)+dzdz,dzdz)
        self.nzf = len(zrange)

        f_mod_fine = np.zeros((self.nzf,self.nt,self.nb))
        for t in range(self.nt):
            for b in range(self.nb):
                spl = InterpolatedUnivariateSpline(self.z_grid, self.f_mod[:,t,b])
                for z in range(self.nzf):
                    f_mod_fine[z,t,b] = spl(zrange[z])

        #flat prior i.e. likelihood not posterior
        like = self.scalefree_flux_likelihood(ff_obs, f_obs_err**2, f_mod_fine)
        #p(z|F)=sum_t[p(F|z,t)p(z,t)]=sum_t[prod_b[N(F_obs_band,F_model_band,sigma_obs_band)]p(z,t)]
        return like.sum()

    def PofZ(self,arr,z,dz):
        """ Returns probability of PZ be at z +-dz/2"""

        self.z_grid = np.load(os.path.join(self.path,'z_grid.npy'))
        self.f_mod = np.load(os.path.join(self.path,'f_mod.npy'))
        (self.nzc,self.nt,self.nb) = np.shape(self.f_mod)

        self.dz = dz
        f_obs,sigma = arr
        self.ngals = len(sigma)

        pofz = np.zeros(self.ngals)
        for n in range(self.ngals):
            pofz[n] = self.getpdf(f_obs[n],sigma[n],z,self.dz)
        return np.array(pofz)

    def cPofZ(self,arr,zx):
        """ Returns cumulative probability of PZ be at z<zx"""
        f_obs,sigma = arr
        #results=np.zeros((self.ngals))#((arr.size))
#         cdfs = np.cumsum(np.array([[self.getpdf(f_obs[n],sigma[n],z,self.dz)
#                         for z in np.arange(0.,zx+self.dz,self.dz)] for n in range(self.ngals)]),
#                         axis=1)
        cdfs = np.cumsum(np.array([self.PofZ(arr,z,self.dz)
                        for z in np.arange(0.,zx+self.dz,self.dz)]),axis=1)
#         mask = pdfs.sum(axis=1)!=0
#         xarrs = arr[mask]['z'][:,np.newaxis] + self.dz
#         masked_pdf = np.where(xarrs<zx, pdfs[mask], 0.)
#         results[mask] = np.sum(masked_pdf, axis=1)
        return cdfs#results

    def getMeanRMS (self,arr):
        """ Returns mean and sqrt variance for
            the photoz pDF, given array
            Note that for assymetric PZ, you are at your
            own risk.
        """
        # in base class we return redshift and zero varinace
        # repeat that here because mean RMS is meaningless for Template SED PDFs
        N=len(arr)
        return arr["z"],np.zeros(N)

    def getMinMax(self,arr):
        """ Returns range of redshifts where p(z) is considerable, i.e.
        no real probability at z<zmin or z>zmax.
        Note that in case of catastrophic outliers, one can have considerable
        amounts of zeros in this range
        """
        # not implemented for Template SED yet
        return arr["z"], arr["z"]

    def NameString(self):
        # not actually sure what this even does
        return "TrueZ"

