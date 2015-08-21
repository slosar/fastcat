#!/usr/bin/env python
import sys

import fastcat
import randomfield as rf
# initialize generator
##fast for debug
#g=fastcat.Generator(grid_spacing_h_Mpc=10.0, smoothing_length_Mpc_h=10.0)
g=fastcat.Generator(grid_spacing_h_Mpc=1.0, smoothing_length_Mpc_h=2.0)
# generate a catalog of objecs
cat=g.genSimple(bias=2,extravar=0.5)
# delete this g to save memory
del g
# Again, but more agressive smoothing required for lognormal
g=fastcat.Generator(grid_spacing_h_Mpc=1.0, smoothing_length_Mpc_h=7.0)
cat2=g.genSimple(bias=2,algorithm="lognormal")
# generate the random catalog
randcat=g.genSimple(algorithm="random")
# make a crappy plot
import pylab
pylab.figure(figsize=(12,12))
pylab.subplot(2,2,1)
pylab.plot (cat["r"], cat["r"]*cat["ra"],'b.',markersize=1.0)
pylab.plot (cat2["r"], cat2["r"]*cat2["ra"],'g.',markersize=1.0)
pylab.plot (randcat["r"], randcat["r"]*randcat["ra"],'r.',markersize=1.0)
pylab.xlabel("z [Mpc/h]")
pylab.ylabel("x [Mpc/h]")
pylab.subplot(2,2,2)
pylab.plot (cat["r"], cat["r"]*cat["dec"],'b.',markersize=1.0)
pylab.plot (cat2["r"], cat2["r"]*cat2["dec"],'g.',markersize=1.0)
pylab.plot (randcat["r"], randcat["r"]*randcat["dec"],'r.',markersize=1.0)
pylab.xlabel("z [Mpc/h]")
pylab.ylabel("y [Mpc/h]")
pylab.subplot(2,2,3)
pylab.plot (cat["ra"], cat["dec"],'b.',markersize=1.0)
pylab.plot (cat2["ra"], cat2["dec"],'g.',markersize=1.0)
pylab.plot (randcat["ra"], randcat["dec"],'r.',markersize=1.0)
pylab.xlabel("ra [rad]")
pylab.ylabel("dec [rad]")

pylab.show()
