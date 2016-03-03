#!/usr/bin/env python
import os
from optparse import OptionParser
from vc_sub import *
import numpy as np 

parser = OptionParser()
if (os.environ.has_key("COLORE_EXEC")):
    cpath=os.environ["COLORE_EXEC"]
else:
    cpath="../CoLoRe/"

opath="/astro/u/anze/Data/colore/"

parser.add_option("--cpath", dest="cpath", default=cpath,
                  help="Path to CoLoRe (will add /CoLoRe for executable)", type="string")
parser.add_option("--outpath", dest="outpath", default=opath,
                  help="Path to output path", type="string")
parser.add_option("--stype", dest="stype", default="wq",
                  help="Submission type (exec,wq)", type="string")
parser.add_option("--N", dest="Nr", default=10,
                  help="Number of realizations", type="int")
parser.add_option("--Ngrid", dest="Ngrid", default=128,
                  help="FFT size", type="int")
parser.add_option("--nodes", dest="parser.add_option("--seed", dest="seed", default=1000,
                  help="Random seed", type="int")

(o, args) = parser.parse_args()

#write distributions 
writeDists(o) 
# now loop over realizations
for i in range(o.Nr):
    execCoLoRe(i,o)


