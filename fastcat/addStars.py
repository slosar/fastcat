#
# adds "stars" to the catalog using a healpix to sample them from
# 
#

import numpy as np
import healpy as hp
import astropy.table as aptable
import catalog as cat
class AddStars:
    """
    Adds objects to a catalog
    """
    default_zmean=0.5
    default_zsigma=0.1
    default_star_file="/project/projectdirs/lsst/LSSWG/Stars/stellar_density_table.fits"

    @staticmethod
    def registerOptions (parser):
        parser.add_option("--stars_zmean", dest="star_zmean", type="float",default=AddStars.default_zmean,
                  help="Mean 'redshift' for stars")
        parser.add_option("--stars_zsigma", dest="star_zsigma", type="float",default=AddStars.default_zsigma,
                  help="Sigma 'redshift' for stars")
        parser.add_option("--stars_catsim_file",dest='star_file',type="string", default=AddStars.default_star_file,
                  help="Input stellar catalog")

    def __init__(self, zmean=default_zmean, zsigma=default_zsigma, star_file=default_star_file, options=None):
        if options is not None:
            self.zmean=options.star_zmean
            self.zsigma=options.star_zsigma
            self.star_file=options.star_file
        else:
            self.zmean=zmean
            self.zsigma=zsigma
            self.star_file=star_file
            
        table=aptable.Table.read(self.star_file)
        self.smap = hp.ma(table['values'].data)
        self.smap.mask = table['mask'].data
        
    def generateStarCatalog(self, nstars,nside2=65536):
        """ 
        Generates  star catalog  by sampling from healpix map.
        nstars : number of stars to generate
        nside2 : healpix grid use to subsample pixels 

        note that since nside2 is never allocated it can be arbitrarily large.
        """

        nside=hp.get_nside(self.smap.data)
        ra=[]
        th=[]
        filled_pixels = np.where(self.smap>0)[0]
        densities = self.smap[filled_pixels]
        kpix = np.random.choice(filled_pixels,size=nstars,p=densities/np.sum(densities))
        bincounts = np.bincount(kpix)
        kpix2 = np.unique(kpix)
        counts=bincounts[bincounts>0]
        hh=nside2**2/nside**2
        i=0
        for i,c in enumerate(counts):
            rpix=np.random.randint(0,high=hh,size=c)
            nestpix=hp.ring2nest(nside,kpix2[i])
            theta, phi = hp.pix2ang(nside2,hh*nestpix+rpix,nest=True)
            theta=90.-theta*180./np.pi
            phi=phi*180./np.pi
            for j in range(0,len(theta)):
                ra.append(phi[j])
                th.append(theta[j])
        ra=np.array(ra)
        dec=np.array(th)
        z=np.random.normal(self.zmean, self.zsigma, nstars)
        catalog=cat.Catalog(nstars)
        catalog['ra']=ra
        catalog['dec']=dec
        catalog['z']=z
        return catalog
            
