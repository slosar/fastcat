#
# Humna window is just a think wrapper around healpix window
# 
#

import numpy as np
import astropy.coordinates as ac
import healpy as hp

from window_Healpix import WindowHealpix

class WindowHumna(WindowHealpix):
    """
    Implements a healpix based map
    """
    typestr='healpix'


    @staticmethod
    def registerOptions (parser):
        parser.add_option("--wf_humnapath", dest="hpath",
                  default="/project/projectdirs/lsst/LSSWG/HumnaDepthVariations",
                  help="Path to humna depth maps", type="string")
        parser.add_option("--wf_humnamap", dest="humnamap", type="string",
                  help="humna map type [nodither, reprandom]", default="nodither")
        parser.add_option("--wf_dlogndmlim", dest="dlogndmlim", type="float",
                  help="change in number of sources as a function of depth",
                  default=0.1)

    def __init__(self, options=None):
        if options is not None:
            if options.humnamap=="nodither":
                mapfn="/coaddM5Data_masked_rBand_NoDither.npz"
            elif options.humnamap=="reprandom":
                mapfn="/coaddM5Data_masked_rBand_RepulsiveRandomDitherFieldPerVisit.npz"
            else:
                print "Unknown humna type map"
                stop()
            mapfn=options.hpath+mapfn
            print "     Reading mapfn",mapfn,"..."
            hmap=np.load(mapfn)
            mask=hmap['mask']
            vals=hmap['metricValues']
            mx=vals.max()
            vals=np.exp(-options.dlogndmlim*(mx-vals))
            vals[mask]=0.0
            amask=np.where(mask==False)
            cmin,cmax,cmean=vals[amask].min(), vals[amask].max(), vals[amask].mean()
            print "Window func min, max, mean:",cmin,cmax,cmean
            info="HumnaDepthVariations map=%s dlogndmlim=%f"%(options.humnamap,options.dlogndmlim)
            shortinfo=options.humnamap+"_"+str(options.dlogndmlim)
            WindowHealpix.__init__(self,vals,info,shortinfo)
        else:
            print "Humna really needs options"
            stop()
            
