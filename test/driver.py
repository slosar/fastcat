#!/usr/bin/env python
import sys
sys.path=["../fastcat"]+sys.path
import fastcat
import astropy.units as u
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--fov", dest="fov", default=3.0,
                  help="Field of view (degrees)", metavar="value", type="float")
parser.add_option("--fast", dest="fast", default=False,
                  action="store_true", help="Settings very fast options for quick test, sets N=10")
parser.add_option("-N", dest="N", default=10000,
                  help="Number of objects to create", metavar="value", type="int")
parser.add_option("--grid_spacing", dest="gspace", default=1,
                  help="Grid Spacing", metavar="value", type="float")
parser.add_option("--smooth", dest="smooth", default=2,
                  help="Smoothing in Mpc/h", metavar="value", type="float")
parser.add_option("--bias", dest="bias", default=2.0,
                  help="bias of tracer", metavar="value", type="float")
parser.add_option("--zmean", dest="zmean", default=1.6,
                  help="Mean redhisft", metavar="value", type="float")
parser.add_option("--deltaz", dest="deltaz", default=0.2,
                  help="variance in redshift", metavar="value", type="float")
parser.add_option("--iesig", dest="iesig", default=0.3,
                  help="intrinsic ellipiticy sigma", metavar="value", type="float")
parser.add_option("--seed", dest="seed", default=123,
                  help="Smoothing in Mpc/h", metavar="value", type="int")
parser.add_option("--algo", dest="alg", default=None,
                  help="Algorithm to use: peaks, lognormal, random", 
                  metavar="value", type="string")
parser.add_option("--phosim", dest="phosim", default=None,
                  help="Set to make a phosim file", 
                  metavar="value", type="string")
parser.add_option("--phosim_header", dest="ps_header", default=None,
                  help="Where to read phosim header from. Empty for no header", 
                  metavar="value", type="string")
parser.add_option("--phosim_many", dest="ps_many", default=False,
                  action="store_true", help="If true, create per obj file")
parser.add_option("--phosim_size", dest="ps_size", default=2.0,
                  help="Size in arcsec of sersic gals", metavar="value", type="float")

(o, args) = parser.parse_args()
if (o.fast):
    o.gspace=10
    o.smooth=10
    o.algo="lognormal"
    o.N=10

maxz=o.zmean+o.deltaz*5
gen=fastcat.Generator(zmax=maxz,size=o.fov*u.deg,grid_spacing_h_Mpc=o.gspace, 
                      smoothing_length_Mpc_h=o.smooth,seed=o.seed)
cat=gen.genSimple(N=o.N,bias=o.bias,zdist=fastcat.ZDist(o.zmean,o.deltaz),
                edist=fastcat.EllipticityDist(o.iesig),algorithm=o.algo)

if o.phosim:
    if o.ps_header:
        header=open(o.ps_header).read()
    else:
        header=""
    cat.dumpPhoSim(o.phosim, header=header, manyFiles=o.ps_many, sedName="../sky/sed_flat.txt", 
                   objtype="sersic2D", ssize=o.ps_size*u.arcsec)
