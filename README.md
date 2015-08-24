# fastcat

Fast and dirty creation of fast and dirty mock galaxy catalogs

This packages takes randomfield to generate a random delta field and
then creates a set of astro objects that have the correct large scale
bias using two different algorithms (or make random catalogs). It can
also produce random ellipticities and should be able to use random
field to create shears at the correct positions. 
For the time being, rmags are drawn from a Gaussian distribution.

## Algorithms

* Peaks: take delta at a field, perturb it with a Gaussian and pick it if it is above
a given limit

* lognormal: make a lognormal transormation of the field and then just sample from it 
in the same way you would from a probability distribution. Need to use much more agressive
smoothing otherwise it gets too concentrated

## Usage

The basic driver is in test/driver.py.

Running with help will tell you things:

```
> test/driver.py --help
Usage: driver.py [options]

Options:
  -h, --help            show this help message and exit
  --fov=value           Field of view (degrees)
  --fast                Settings very fast options for quick test, sets N=10
  -N value              Number of objects to create
  --grid_spacing=value  Grid Spacing
  --smooth=value        Smoothing in Mpc/h
  --bias=value          bias of tracer
  --zmean=value         Mean redhisft
  --deltaz=value        variance in redshift
  --iesig=value         intrinsic ellipiticy sigma
  --seed=value          Smoothing in Mpc/h
  --algo=value          Algorithm to use: peaks, lognormal, random
  --phosim=value        Set to make a phosim file
  --phosim_header=value
                        Where to read phosim header from. Empty for no header
  --phosim_many         If true, create per obj file
  --phosim_size=value   Size in arcsec of sersic gals
  --h5read=value        Instead of creating dataset, read it from H5
  --h5write=value       Write to hdf5 file specified on command line
```

For more information, best look at https://github.com/DarkEnergyScienceCollaboration/BremertonRoundTrip/

For direct usage, it is all very easy, see example in test/test.py

The basic procedure is the following:

* Set up generator object. This will use randomfields to generate delta field

* Use genSimple to generate catalog

* Catalog is is of fastcat.Catalog object as written out in fastcat/catalog.py. 
  It holds a structured array which you can access directly.
  Eg. cat["ra"] will give you 1D array of ra coordinas. Valid names are:

  * "ra", "dec": float, ra,dec coordinates in radians
  * "z": float, redshift
  * "r": float, distance in Mpc/h
  * "rmag": float, rmagnitude
  * "e1","e2" : float, intrinsic ellipticity
  * "g1","g2" : float, shears

## TODO

* Shears not yet implemented.



