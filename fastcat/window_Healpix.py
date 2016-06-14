import numpy as np
import astropy.coordinates as ac
import healpy as hp

from window_Base import WindowBase

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
        dset.attrs['short_info']=self.sinfo
        dset.attrs['info']=self.info
        

    @staticmethod
    def readH5 (dataset):
        """ Tries to read from H5.
            If not matched, return None
        """
        if dataset.attrs['type']==WindowHealpix.typestr:
            info=dataset.attrs['info']
            try:
                sinfo=dataset.attrs['short_info']
            except KeyError:
                sinfo=info
            hmap=dataset.value
            return WindowHealpix(hmap,info,sinfo)
        else:
            return None
        
    def NameString (self):
        return "hpix_"+self.sinfo
        
