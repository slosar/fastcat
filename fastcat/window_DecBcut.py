import numpy as np
import astropy.coordinates as ac

from window_Base import WindowBase

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

