from __future__ import print_function, division

try:
    import treecorr
except:
    pass

import math

class TCDriver(object):
    """ 
    Provides convenient interface to TreeCorr correlation function package.
   """
    
    def __init__ (self, catalog, randcatalog, min_sep=1, Nbins=100):

        try:
            dir(treecorr)
        except NameError:
            print ("Don't have tree corr.")
            stop()
        self.catalog=catalog
        self.randcatalog=randcatalog
        self.max_sep=self.catalog.max_sep*60 ## Let's work with arcmins here
        self.min_sep=min_sep
        self.bin_size=math.log(self.max_sep/self.min_sep)/Nbins



    def NN3DCorrelation(self, min_sep=1, max_sep=200, bin_size=0.5):
        """
        Caclulates 3D correlation function of objects using Catalog's ra, dec, r
        Requires randcatalog to exist. Distance units are Mpc/h.

        Returns tuple (logr, meanlogr, xi, xivar)
        
        """
        catS = treecorr.Catalog(ra=self.catalog["ra"], dec=self.catalog["dec"], 
                                r=self.catalog["r"],
                                ra_units="radians", dec_units="radians")
        if (self.randcatalog):
            catR = treecorr.Catalog(ra=self.randcatalog["ra"], 
                                    dec=self.randcatalog["dec"], 
                                    r=self.randcatalog["r"], ra_units="radians", 
                                    dec_units="radians")
        else:
            print ("Need random catalog for NN")
            stop()
        dd=treecorr.NNCorrelation(min_sep=min_sep, bin_size=bin_size, 
                                  max_sep=max_sep)
        dr=treecorr.NNCorrelation(min_sep=min_sep, bin_size=bin_size, 
                                  max_sep=max_sep)
        rr=treecorr.NNCorrelation(min_sep=min_sep, bin_size=bin_size, 
                                  max_sep=max_sep)
            
        dd.process(catS)
        dr.process(catS, catR)
        rr.process(catR)
        xi, xivar=dd.calculateXi(rr,dr)
        logr = dd.logr 
        meanlogr= dd.meanlogr
        return (logr, meanlogr, xi, xivar)


    def NNCorrelation(self):
        """
        Caclulates 2D correlation function using Catalog's ra, dec.
        Requires randcatalog to exist.

        Returns tuple (logr, meanlogr, xi, xivar)
        
        """
        catS = treecorr.Catalog(ra=self.catalog["ra"], dec=self.catalog["dec"], 
                                ra_units="radians", dec_units="radians")
        if (self.randcatalog):
            catR = treecorr.Catalog(ra=self.randcatalog["ra"], 
                                    dec=self.randcatalog["dec"], ra_units="radians", 
                                    dec_units="radians")
        else:
            print ("Need random catalog for NN")
            stop()
        dd=treecorr.NNCorrelation(min_sep=self.min_sep, bin_size=self.bin_size, 
                                  max_sep=self.max_sep, sep_units='arcmin')
        dr=treecorr.NNCorrelation(min_sep=self.min_sep, bin_size=self.bin_size,
                                  max_sep=self.max_sep, sep_units='arcmin' )
        rr=treecorr.NNCorrelation(min_sep=self.min_sep, bin_size=self.bin_size,
                                  max_sep=self.max_sep, sep_units='arcmin')
        dd.process(catS)
        dr.process(catS, catR)
        rr.process(catR)
        xi, xivar=dd.calculateXi(rr,dr)
        logr = dd.logr 
        meanlogr= dd.meanlogr
        return (logr, meanlogr, xi, xivar)

    def GGCorrelation(self):
        """
        Caclulates 2D correlation function using Catalog's ra, dec.
        Requires randcatalog to exist.

        Returns tuple (logr, meanlogr, xip, xim, xivar)
        """
        catS = treecorr.Catalog(ra=self.catalog["ra"], dec=self.catalog["dec"], 
                  ra_units="radians", dec_units="radians", g1=self.catalog["g1"],
                  g2=self.catalog["g2"] )
        dd=treecorr.GGCorrelation(min_sep=self.min_sep, bin_size=self.bin_size, 
                                  max_sep=self.max_sep, sep_units='arcmin')

        dd.process(catS)
        logr = dd.logr 
        meanlogr = dd.logr
        xip=dd.xip
        xim=dd.xim
        xivar=dd.varxi
        return (logr, meanlogr, xip, xim, xivar)

    def NGCorrelation(self):
        """
        Caclulates 2D correlation function using Catalog's ra, dec.
        Requires randcatalog to exist.

        Returns tuple (logr, meanlogr, xi, xi_im, xivar)
        """
        catN = treecorr.Catalog(ra=self.catalog["ra"], dec=self.catalog["dec"], 
                  ra_units="radians", dec_units="radians", g1=self.catalog["g1"],
                  g2=self.catalog["g2"] )
        dd=treecorr.NGCorrelation(min_sep=self.min_sep, bin_size=self.bin_size, 
                                  max_sep=self.max_sep, sep_units='arcmin')

        dd.process(catN,catN)
        logr = dd.logr 
        meanlogr = dd.logr
        xi=dd.xi
        xi_im=dd.xi_im
        xivar=dd.varxi
        return (logr, meanlogr, xi, xi_im, xivar)


                                
