from __future__ import print_function, division

import numpy as np

class Catalog(object):
    """ Basic object to hold a catalog of observed astronomical objects.
        Intentially very simple for the time being"""
    
    def __init__ (self, N,meta=None):
        self.data=np.zeros(N,dtype=[('ra',float),('dec',float),('z',float),('r',float),
                                     ('rmag',float), ('e1',float), ('e2',float),
                                     ('g1',float),('g2',float)])
        self.meta=meta

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key,item):
        self.data[key]=item
