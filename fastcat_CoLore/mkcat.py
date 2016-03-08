#!/usr/bin/env python
import sys
from subs import *
sys.path=["fastcat"]+sys.path
import catalog as cat
import numpy as np 
from optparse import OptionParser
import astropy.coordinates as ac
import datetime

version='0.1'

ipath="/astro/u/anze/Data/colore/"
opath="/astro/u/anze/Data/colcat/"
parser = OptionParser()
parser.add_option("--ipath", dest="ipath", default=ipath,
                  help="Path to colore output", type="string")
parser.add_option("--opath", dest="opath", default=opath,
                  help="Path to colore output", type="string")
parser.add_option("--N", dest="Nr", default=10,
                  help="Number of realizations", type="int")
parser.add_option("--decmin", dest="decmin", default=-70,
                  help="Number of realizations", type="float")
parser.add_option("--decmax", dest="decmax", default=10,
                  help="Number of realizations", type="float")
parser.add_option("--bcut", dest="bcut", default=5,
                  help="Cut in galactic b", type="float")
(o, args) = parser.parse_args()

bz=np.genfromtxt(o.ipath+'/bz.txt', dtype=None, names=["z","bz"])
dNdz=np.genfromtxt(o.ipath+'/Nz.txt', dtype=None, names=["z","dNdz"])

for i in range(o.Nr):
    print "Reading set ",i
    gals,inif=readColore(o.ipath+"/Set%i"%i)
    print len(gals)," galaxies read."
    if (len(gals)==0):
        print "No galaxies!"
        stop()
    gals=gals[(gals['DEC']>o.decmin) & (gals['DEC']<o.decmax)]
    print len(gals)," after dec cut"
    if hasattr(ac,"SkyCoord"):
        galb=np.array(ac.SkyCoord(ra=gals['RA'], dec=gals['DEC'], frame='icrs', unit='deg').galactic.b)
    else:
        print "You have an old astropy! This will be slow!!"
        #old school
        cra=gals["RA"]
        cdec=gals["DEC"]
        galb=np.zeros(len(gals))
        for i in xrange(len(gals)):
            galb[i]=ac.ICRSCoordinates(ra=cra[i], dec=cdec[i], unit=('deg','deg')).galactic.b.degrees

    gals=gals[abs(galb)>o.bcut]
    print len(gals)," after b cut"
    N=len(gals)
    meta={}
    for k,v in inif.items():
        meta['colore_'+k]=v
    meta['version']=version
    meta['timestamp']=str(datetime.datetime.now())
    window={"type":"decbcut",
            "dec_min":o.decmin,
            "dec_max":o.decmax,
            "b_cut":o.bcut}
    photoz={"type":"fixed_sigma",
            "sigma":"0.01"}
    
    cat=cat.Catalog(N, fields=['ra','dec','z_real_t','z_red_t','z_error'],dNdz=dNdz, bz=bz,
                    photoz=photoz,window=window,meta=meta)
    cat['ra']=gals['RA']
    cat['dec']=gals['DEC']
    cat['z_real_t']=gals['Z_COSMO']
    cat['z_red_t']=gals['Z_COSMO']+gals['DZ_RSD']
    cat['z_error']=np.random.normal(0,photoz["sigma"],N)
    cat.dumpH5(o.opath+'/catalog%i.h5'%(i))


