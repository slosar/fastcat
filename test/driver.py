#!/usr/bin/env python
import sys
sys.path=["../fastcat"]+sys.path
import fastcat
import astropy.units as u
from optparse import OptionParser
import math

parser = OptionParser()
parser.add_option("--fov", dest="fov", default=3.0,
                  help="Field of view (degrees)", metavar="value", type="float")
parser.add_option("--fast", dest="fast", default=False,
                  action="store_true", help="Settings very fast options for quick test, sets N=10")
parser.add_option("--fast2", dest="fast2", default=False,
                  action="store_true", help="Settings very medium speed, sets N=10000")
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
parser.add_option("--deltaz", dest="deltaz", default=0.05,
                  help="variance in redshift", metavar="value", type="float")
parser.add_option("--iesig", dest="iesig", default=0.3,
                  help="intrinsic ellipiticy sigma", metavar="value", type="float")
parser.add_option("--seed", dest="seed", default=123,
                  help="Smoothing in Mpc/h", metavar="value", type="int")
parser.add_option("--boxpad", dest="boxpad", default=1.5,
                  help="What factor to pad the box with", metavar="value", type="float")
parser.add_option("--algo", dest="algo", default="peaks",
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
parser.add_option("--h5read", dest="h5read", default=None,
                  help="Instead of creating dataset, read it from H5", 
                  metavar="value", type="string")
parser.add_option("--h5write", dest="h5write", default=None,
                  help="Write to hdf5 file specified on command line", 
                  metavar="value", type="string")
parser.add_option("--treecorr", dest="treecorr", default=None,
                  help="Set to use TreeCorr to calculate corr functions. "
                  " Output to filename specified")
parser.add_option("--treecorr3D", dest="treecorr3D", default=None,
                  help="Set to use TreeCorr to calculate 3D corr function of "
                  " tracers. Output to filename specified")

(o, args) = parser.parse_args()
if (o.fast):
    o.gspace=10
    o.smooth=10
    o.algo="lognormal"
    o.N=10

if (o.fast2):
    o.gspace=3
    o.smooth=6
    o.N=10000

if (o.h5read):
    cat=fastcat.Catalog(0)
    cat.readH5(o.h5read)
else:
    maxz=o.zmean+o.deltaz*5
    gen=fastcat.Generator(zmax=maxz,size=o.fov*u.deg,grid_spacing_h_Mpc=o.gspace, 
                boxpad=o.boxpad, smoothing_length_Mpc_h=o.smooth,seed=o.seed)
    cat=gen.genSimple(N=o.N,bias=o.bias,zdist=fastcat.ZDist(o.zmean,o.deltaz),
                edist=fastcat.EllipticityDist(o.iesig),algorithm=o.algo)


if (o.treecorr):
    ## first need to generate random catalog
    catr=gen.genSimple(N=2*o.N,bias=o.bias,zdist=fastcat.ZDist(o.zmean,o.deltaz),
                       edist=fastcat.EllipticityDist(o.iesig),algorithm="random")
    ## now calculate 
    tcd=fastcat.TCDriver(cat,catr)
    logr, meanlogr, xinn, xivarnn = tcd.NNCorrelation()
    logr, meanlogr, xip, xim, xivargg = tcd.GGCorrelation()
    logr, meanlogr, xing, xingi, xivarng = tcd.NGCorrelation()
    of=open(o.treecorr,'w')
    of.write('# logr_nominal logr xi_dd xi_gg+ xi_gg- xi_dg xi_dg_imag  (each xi is value error)\n')
    for  a,b,m1,e1, m2,m3,e23, m4, m5, e45 in zip(logr, meanlogr, xinn, xivarnn,
                                                  xip, xim, xivargg, xing, xingi, xivarng):
        of.write("%g %g %g %g %g %g %g %g %g %g %g %g \n "%(math.exp(a),math.exp(b),m1,e1,
                                                            m2,e23,m3,e23,m4,e45,m5,e45))
    of.close()


if (o.treecorr3D):
    ## first need to generate random catalog
    catr=gen.genSimple(N=2*o.N,bias=o.bias,zdist=fastcat.ZDist(o.zmean,o.deltaz),
                       edist=fastcat.EllipticityDist(o.iesig),algorithm="random")
    ## now calculate 
    tcd=fastcat.TCDriver(cat,catr)
    logr, meanlogr, xinn, xivarnn = tcd.NN3DCorrelation(1,1000,0.05)
    of=open(o.treecorr3D,'w')
    of.write('# logr_nominal logr xi_dd xi_dd_error\n')
    for  a,b,m1,e1 in zip(logr, meanlogr, xinn, xivarnn):
        of.write("%g %g %g %g \n "%(math.exp(a),math.exp(b),m1,e1))
    of.close()


if o.phosim:
    if o.ps_header:
        header=open(o.ps_header).read()
    else:
        header=""
    cat.dumpPhoSim(o.phosim, header=header, manyFiles=o.ps_many, sedName="../sky/sed_flat.txt", 
                   objtype="sersic2D", ssize=o.ps_size*u.arcsec)
if o.h5write:
    cat.dumpH5(o.h5write)

