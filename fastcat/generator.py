from __future__ import print_function, division

import sys, math, random
import numpy as np
import astropy.units as u
import randomfield as rf
import randomfield.cosmotools as ct
import randomfield.powertools as pt
import catalog as cat
from scipy.interpolate import interp1d
from scipy.optimize import brentq
from numpy import sqrt, exp


class ZDist(object):
    """ Characterises the redshift distribution with a Gaussian. """
    def __init__ (self,zmean=1.6, deltaz=0.2):
        self.zmean=zmean
        self.deltaz=deltaz

    def pofz(self,z):
        return exp(-(z-self.zmean)**2/(2*self.deltaz**2))

    def refz(self):
        return self.zmean

    def sampleZ(self):
        return np.random.normal(self.zmean,self.deltaz)

class EllipticityDist(object):
    """ Characterises the redshift distribution with a Gaussian. """
    def __init__ (self, abs_ellipticity_sigma=0.3):
        self.eonesigma=abs_ellipticity_sigma**2/2 ## divide by two per component
        
    def sampleEOne(self):
        return np.random.normal(0.0,self.eonesigma)


class MagDist(object):
    """ Characterises the redshift distribution with a Gaussian. """
    def __init__ (self, rmag=19, deltarmag=2):
        self.rmagmean=19
        self.deltarmag=deltarmag
        
    def sampleRMag(self,bias):
        ## in principle, more biased object are
        ## brigher, so we can simplysim this. But so far
        ## we ignore it
        return np.random.normal(self.rmagmean,self.deltarmag)


class Generator(object):
    def __init__ (self,zmax=2.5, size=2.0*u.deg, grid_spacing_h_Mpc=1.0, smoothing_length_Mpc_h=2.0,
                  seed=123,boxpad=1.2,cosmology=None, power=None):
        if cosmology is None:
            self.cosmology = ct.create_cosmology()
            if power is None:
                power = pt.load_default_power()
        else:
            self.cosmology = cosmology
            if power is None:
                power = ct.calculate_power(
                    self.cosmology, k_min=self.k_min, k_max=self.k_max,
                    scaled_by_h=True)
        self.power = pt.validate_power(power)
        ## first work out how big a box do we need
        dmax=self.cosmology.comoving_distance(zmax)*self.cosmology.h/u.Mpc
        self.Nz=int(dmax/grid_spacing_h_Mpc*boxpad)+1
        self.thetao2=(size/u.rad)/2.0
        self.Nx=int(dmax*self.thetao2*2/grid_spacing_h_Mpc*boxpad)+1

        ## actually we need things to be divisible by 4
        while (self.Nz%4):
            self.Nz+=1
        while (self.Nx%4):
            self.Nx+=1
        self.Ny=self.Nx

        memsize=(self.Nx*self.Ny*self.Nz*8/1e9)
        print("Will need box of size ",self.Nx,"x",self.Ny,"x",self.Nz, "; about ",memsize,"Gb")
        if (memsize>4.):
            print("If this is too much for your computer, ask your mom for a bigger one.")

        self.gen=rf.Generator(self.Nx, self.Ny, self.Nz, grid_spacing_Mpc_h=grid_spacing_h_Mpc,
                                       cosmology=cosmology, power=power, verbose=True)

        self.delta=self.gen.generate_delta_field(
            smoothing_length_Mpc_h=smoothing_length_Mpc_h, seed=seed)
        self.var=self.delta.var()
        self.zs=ct.get_redshifts(self.cosmology, self.delta,grid_spacing_h_Mpc)
        ## get interpolator from index to z
        self.k2z=interp1d(range(self.Nz),self.zs[0,0,:])
        self.z2k=interp1d(self.zs[0,0,:],range(self.Nz))

        self.grid_spacing_h_Mpc=grid_spacing_h_Mpc
        
    def genSimple(self, N=10000, bias=2.0, zdist=ZDist(), rmagdist=MagDist(),
                  edist=EllipticityDist(), algorithm="peaks", extravar=0.5):

        """
        Generate a catalog of fake objects

        Parameters
        ----------
        N : int
            Number of objects to generata
        bias : float
            The bias at reference redshift (as specified by zdist.refz()) 
            The algorithm actually operates with z=0 delta field so we convert
            this bias factor to bias factor today.
        zdist: object
            Object specifying redshift window function (Gaussian for the time being)
        rmagdist: object
            Object specifying rmag distribution (Gaussian for the time being)
        edist: object
            Object specifying intrinsic ellipticity distribution
        algorithm: string
            How to distribute objects. Valid options are i) "peaks" (put galaxy if
            delta value + extra white noise > limit determined so that bias is right).
            ii) "lognormal" probabilistically sample from lognormal transformation 
            (with the right bias), iii) "random" for random catalogs (just window function)

        extravar: float
            Extra variance if algorithm=="peaks"
        """

        random=False
        peaks=True
        if algorithm=="peaks":
            pass
        elif algorithm=="lognormal":
            peaks=False
        elif algorithm=="random":
            random=True
        else:
            print ("Bad algorithm")
            stop()

        ## we have deltafield at z=0, so first get bias right
        # now convert bias at effective z to bias now
        gf=ct.get_growth_function(self.cosmology,np.array([0,zdist.refz()]))
        biasnow=bias*gf[1]/gf[0]
        
        if (peaks):
            cg0=(1.+extravar)
            #Now need to calculate clipping at which I will produce sufficient number of objects
            # from maple
            fun= lambda l: (-(math.exp(-l ** 2 / cg0 ** 2 / 0.2e1) * math.sqrt(0.2e1) * 
                            math.pi ** (-0.1e1 / 0.2e1) / cg0 / (math.erf(math.sqrt(0.2e1) 
                           / cg0 * l / 0.2e1) - 0.1e1))-biasnow)
            ## first healtcheck
            if ((fun(-3)>0) or (fun(3)<0)):
                print ("You are asking for a bias which is a bit extreme.")
                print (fun(-3)+biasnow, fun(3)+biasnow)
                stop()
            lim=brentq(fun,-3.0, 3.0)
            ## lim is now cut in sources of deltafield + extra*var. Let us just go through all of them
            lim*=sqrt(self.var)
            esg=sqrt(self.var*extravar)
            print ("Picking up above ",lim, " with scatter ",esg," sig=",sqrt(self.var))

        else:
            lognmax=exp(1+self.delta.max()*biasnow)
            print ("lognmax=",lognmax)
        tx=self.Nx/2
        ty=self.Ny/2

        zlist=[]
        dlist=[]
        ralist=[]
        rlist=[]
        klist=[]
        thetao2sq=self.thetao2**2
        tcc=0
        for cc in xrange(N):
            while(True):
                tcc+=1
                z=zdist.sampleZ()
                kre=self.z2k(z) ## radial distance in indices
                ## now sample radially
                r=sqrt(np.random.uniform()*thetao2sq)
                phi=np.random.uniform(0,2*math.pi);
                ra=r*np.sin(phi)
                dec=r*np.cos(phi)
                
                k=int(kre)
                i=int(ra*kre+tx)
                j=int(dec*kre+ty)
                accept=random
                if (not accept):
                    delta=self.delta[i,j,k]
                    if (peaks):
                        if esg>0:
                            accept=(delta+np.random.normal(0,esg))>lim
                        else:
                            accept=delta>lim
                    else:
                        accept=(exp(1+delta*biasnow)/lognmax)>np.random.uniform()

                ## now try to sample
                if accept:
                    ## ok, now we really have it.
                    zlist.append(z)
                    dlist.append(dec)
                    ralist.append(ra)
                    rlist.append(kre*self.grid_spacing_h_Mpc)
                    break

        print ("Proposal accept rate:",float(N)/tcc)
        toret=cat.Catalog(N)
        toret["ra"]=np.array(ralist)
        toret["dec"]=np.array(dlist)
        toret["z"]=np.array(zlist)
        toret["r"]=np.array(rlist)*self.grid_spacing_h_Mpc
        toret["rmag"]=np.array([rmagdist.sampleRMag(bias) for i in xrange(N)])
        toret["e1"]=np.array([edist.sampleEOne() for i in xrange(N)])
        toret["e2"]=np.array([edist.sampleEOne() for i in xrange(N)])
        toret["g1"]=np.zeros(N)
        toret["g2"]=np.zeros(N)
        
        return toret

    

