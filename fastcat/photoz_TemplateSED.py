#
# Template SED photo-z object
#
import numpy as np
from numpy.lib import recfunctions
import os

import sys
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d

from photoz_Base import PhotoZBase

# filterloc = 'Users/aimalz/Photoz-tools/filters'
# sedloc = 'Users/aimalz/Photoz-tools/seds'

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
        #store from Boris' notebook and read in here:
        #flux noise sigma_obs (fix this as a constant),
        #flux-z model F_model (3D array, write out from Boris' notebook, then read in and interpolate),
        #p(z,t) (or just give likelihood)

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

        #what specifies a PDF in Template SED method?
        #the fluxes and noise corresponding to the bands
        #which means they need to be simulated here

        ztrue = arr['z']

        #select a template
        templates = ['El_B2004a.sed']+['Sbc_B2004a.sed','Scd_B2004a.sed']
        templates = templates +['Im_B2004a.sed','SB3_B2004a.sed','SB2_B2004a.sed','ssp_25Myr_z008.sed','ssp_5Myr_z008.sed']

        #read in f_mod files, interpolate, get values of f_mod_b
        self.z_grid = np.load(os.path.join(self.path,'z_grid.npy'))
        self.f_mod = np.load(os.path.join(self.path,'f_mod.npy'))
        (nz,nt,nb) = np.shape(self.f_mod)
        self.nt = nt
        self.nb = nb
        #self.nz = nz
        print('nz='+str(nz)+', z_grid='+str(len(self.z_grid)))

        ngals = len(ztrue)
        f_mod_o = np.zeros((self.nb, ngals))
        for z in range(ngals):
            templateno = np.random.choice(range(self.nt))
            for b in range(self.nb):#check this manipulation
                spl = InterpolatedUnivariateSpline(self.z_grid, self.f_mod[:,templateno,b])
                f_mod_o[b][z] = spl(ztrue[z])
        print('f_mod_o'+str(np.shape(f_mod_o)))

        #select sigma_b - 1% for now
        sigma = 0.01*f_mod_o

        #select observed fluxes f_obs_b = f_mod_b + sigma_b*rando
        f_obs = np.zeros((self.nb, ngals))#np.zeros(self.nb)
        for n in xrange(len(ztrue)):
            for b in range(self.nb):
                f_obs[b][n] = f_mod_o[b][n] + sigma[b][n] * np.random.normal()
        print('f_obs'+str(np.shape(f_obs)))

        #eventually want this to be meaningful output
        return f_obs,sigma

    def scalefree_flux_likelihood(self, f_obs, f_obs_var, f_mod):
        f_obs = f_obs/np.mean(f_obs)
        f_obs_var = f_obs_var/np.mean(f_obs_var)
        f_mod = f_mod/np.mean(f_mod)

#         lf_obs = np.log(f_obs)
#         lf_obs_var = np.log(f_obs_var)
#         lf_mod = np.log(f_mod)

        var = f_obs_var  # nz * nt * nf
        invvar = np.where(f_obs/var < 1e-6, 0.0, var**-1.0)  # nz * nt * nb
        FOT = np.sum(f_mod * f_obs * invvar, axis=2)  # nz * nt
        FTT = np.sum(f_mod**2 * invvar, axis=2)  # nz * nt
        FOO = np.dot(invvar, f_obs**2)  # nz * nt
        chi2 = FOO - FOT**2.0 / FTT  # nz * nt
        like = np.exp(-0.5*chi2) / np.sqrt(FTT)  # nz * nt
        return like

    def getpdf(self,f_obs,sigma,z,dz):

        ff_obs = (10.**(-.4*f_obs))
        #print('ff_obs '+str(np.shape(ff_obs)))
        f_obs_err = (10.**(.4*np.abs(sigma))-1.)*ff_obs
        #print('f_obs_err '+str(np.shape(f_obs_err)))

        #p(z|F)=sum_t[p(F|z,t)p(z,t)]=sum_t[prod_b[N(F_obs_band,F_model_band,sigma_obs_band)]p(z,t)]
#         summed_t = 0.
#         for t in range(self.nt):
#             prod_b = 1.
#             for b in range(self.nb):
#                 summed_z = 0.
#                 spl = InterpolatedUnivariateSpline(self.z_mod, self.f_mod[t][b])
#                 for z in range(self.nz):
#                     summed_z += spl(zrange[z])
#                 prod_b *= summed_z
#             summed_t += prod_b

        self.z_grid = np.load(os.path.join(self.path,'z_grid.npy'))
        self.f_mod = np.load(os.path.join(self.path,'f_mod.npy'))
        (nz,nt,nb) = np.shape(self.f_mod)
        self.nt = nt
        self.nb = nb
        #self.nz = nz

        zrange = np.arange(max(0.,z-dz/2.),min(z+dz/2.,1.4)+dz**2,dz**2)
        self.nz = len(zrange)
        #print('why me?'+str(self.nz))

        f_mod_fine = np.zeros((self.nz,self.nt,self.nb))
        for t in range(self.nt):
            for b in range(self.nb):
                spl = InterpolatedUnivariateSpline(self.z_grid, self.f_mod[:,t,b])
                for z in range(self.nz):
                    f_mod_fine[z,t,b] = spl(zrange[z])

        #flat prior i.e. likelihood not posterior
        like = self.scalefree_flux_likelihood(ff_obs, f_obs_err**2, f_mod_fine)
        #sum_t[prod_b[N(F_obs_band,F_model_band,sigma_obs_band)]]
        loglike = np.log(like)
        loglike_t = loglike.sum(axis=0)
        like_t = np.exp(loglike_t)
        out = np.sum(like_t)
        return out

    def PofZ(self,arr,z,dz):
        """ Returns probability of PZ be at z +-dz/2"""

#         results=np.zeros((arr.size))
#         pdfs = self.getpdf(z,dz)
#         mask = pdfs.sum(axis=1)!=0
#         xarrs = arr[mask]['z'][:,np.newaxis] + self.dz
#         masked_pdf = np.where(np.abs(xarrs-z)<dz/2., pdfs[mask], 0.)
#         results[mask] = np.sum(masked_pdf, axis=1)
        #N = len(arr)
        #P = np.zeros(N)
        #dzo2 = dz/2
        #P[np.where(abs(arr["z"]-z)<dzo2)] = self.getpdf(z,dz)
        f_obs,sigma = arr
        ngals = len(sigma)

        pofz = []
        for n in range(ngals):
            pofz.append(self.getpdf(f_obs[n],sigma[n],z,dz))

        #print('pofz gives '+str(np.shape(np.array(pofz))))
        return np.array(pofz)#P#results

    def cPofZ(self,arr,zx):
        """ Returns cumulative probability of PZ be at z<zx"""
        results=np.zeros((arr.size))
        pdfs = self.getpdf(arr)
        mask = pdfs.sum(axis=1)!=0
        xarrs = arr[mask]['z'][:,np.newaxis] + self.dz
        masked_pdf = np.where(xarrs<zx, pdfs[mask], 0.)
        results[mask] = np.sum(masked_pdf, axis=1)
        return results

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

    def NameString(self):
        return "TrueZ"

