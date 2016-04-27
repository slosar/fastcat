##
## Class describing window function
##
import numpy as np
import astropy.coordinates as ac
import healpy as hp

## this method does not need an object    
def readWindowH5(dataset):
    name=dataset.attrs['type']
    ## loop over possible types
    ## (why can't i do 'for t in [WindowBase, WindowDecBcut]:'?)
    toret=WindowBase.readH5(dataset)
    if toret is not None: return toret
    toret=WindowDecBcut.readH5(dataset)
    if toret is not None: return toret
    toret=WindowHealpix.readH5(dataset)
    if toret is not None: return toret
    print "Unknown window type!"
    stop()


class WindowBase(object):
    """
    Basic object for describing window function of survey.
    It needs to implement one function: 

    float __call__(ra,dec),

    which returns relative angular probability of seeing
    an object (between 0 and 1).

    ras and decs can also be numpy arrays in which case
    it returns a numpyu array
    """
    typestr="base"
    def __init__(self):
        pass
        
    def __call__(self,ra,dec):
        if (type(ra)==type(1.0)):
            return 1.0
        else:
            return np.ones(len(ra))

    def writeH5 (self,of):
        dset=of.create_dataset("window",data=[])
        dset.attrs['type']=self.typestr

    @staticmethod
    def readH5 (dataset):
        """ Tries to read from H5.
            If not matched, return None
        """
        if dataset.attrs['type']==WindowBase.typestr:
            return WindowBase()
        else:
            return None

    def applyWindow (self,arr):
        """ Applies window function to
            a catalog's array by proportionally (randomly) cutting
            entries given the window efficiency.
        """
        ra=arr['ra']
        dec=arr['dec']
        N=len(ra)
        wf=self(ra,dec)
        mask=(wf>np.random.uniform(0,1,N))
        arr=arr[mask]
        print "After window: ",N,"->",len(arr)
        return arr

    def NameString (self):
        return "FullSky"
    
class WindowDecBcut(WindowBase):
    """
    Implements a basic declination cut and galactic b cut
    """
    typestr='decbcut'
    
    def __init__(self,dec_min,dec_max, b_cut):
        self.dec_min=dec_min
        self.dec_max=dec_max
        self.b_cut=b_cut

    def __call__(self,ra,dec):
        
        if (type(ra)==type(1.0)):
            if (dec<dec_min): return 0.0
            if (dec>dec_max): return 0.0
            galb=np.array(ac.SkyCoord(ra=ra, dec=dec, frame='icrs', unit='deg').galactic.b)
            if (abs(galb)<self.b_cut): return 0.0
            return 1.0
        else:
            toret=np.ones(len(ra))
            toret[np.where(dec<self.dec_min)]=0.0
            toret[np.where(dec>self.dec_max)]=0.0
            galb=np.array(ac.SkyCoord(ra=ra, dec=dec, frame='icrs', unit='deg').galactic.b)
            toret[np.where(abs(galb)<self.b_cut)]=0.0
            return toret

    def writeH5 (self,of):
        dataset=of.create_dataset("window",data=[])
        dataset.attrs['type']=self.typestr
        dataset.attrs['dec_min']=self.dec_min
        dataset.attrs['dec_max']=self.dec_max
        dataset.attrs['b_cut']=self.b_cut

    @staticmethod
    def readH5 (dataset):
        """ Tries to read from H5.
            If not matched, return None
        """
        if dataset.attrs['type']==WindowDecBcut.typestr:
            dec_min=float(dataset.attrs['dec_min'])
            dec_max=float(dataset.attrs['dec_max'])
            b_cut=float(dataset.attrs['b_cut'])
            return WindowDecBcut(dec_min, dec_max,b_cut)
        else:
            return None
        
    def NameString (self):
        return "RaDecCut"


class WindowHealpix(WindowBase):
    """
    Implements a healpix based map
    """
    typestr='healpix'
    
    def __init__(self,healpixmap, info, shortinfo):
        self.map=healpixmap
        self.info=info
        self.sinfo=shortinfo
        self.nside=hp.pixelfunc.npix2nside(len(self.map))
        
    def __call__(self,ra,dec):
        phi=ra/180.*np.pi
        theta=np.pi/2-dec/180.*np.pi
        ndx=hp.pixelfunc.ang2pix(self.nside, theta,phi)
        return self.map[ndx]

    def writeH5 (self,of):
        dset=of.create_dataset("window",data=self.map)
        dset.attrs['type']=self.typestr
        dset.attrs['info']=self.info

    @staticmethod
    def readH5 (dataset):
        """ Tries to read from H5.
            If not matched, return None
        """
        if dataset.attrs['type']==WindowHealpix.typestr:
            info=dataset.attrs['info']
            hmap=dataset.value
            return WindowHealpix(hmap,info)
        else:
            return None
        
    def NameString (self):
        return "hpix_"+self.sinfo
        
