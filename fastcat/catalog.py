from __future__ import print_function, division

import numpy as np

class Catalog(object):
    """ 
    Basic object to hold a catalog of observed astronomical objects.
    Intentially very simple for the time being.

    It holds a structured array which you can access directly.
    Eg. cat["ra"] will give you 1D array of ra coordinas. Valid names are:
    
    * "ra", "dec": float, ra,dec coordinates in radians
    * "z": float, redshift
    * "r": float, distance in Mpc/h
    * "rmag": float, rmagnitude
    * "e1","e2" : float, intrinsic ellipticity
    * "g1","g2" : float, shears
    It also has a placeholder for meta-data, which is empty at the moment
   """
    
    def __init__ (self, N,meta=None):
        self.data=np.zeros(N,dtype=[('ra',float),('dec',float),('z',float),('r',float),
                                     ('rmag',float), ('e1',float), ('e2',float),
                                     ('g1',float),('g2',float)])
        self.meta=meta

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key,item):
        self.data[key]=item

    def dumpPhoSim(self, fname):
        """
        Writes out catalog in format that phosim can chew.

        Parameters
        ----------
        fname : string
        Filename to dump it to. 
        """
        print ("Implement dump phosim.")
        stop()
        
