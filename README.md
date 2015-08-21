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

It is all very easy, see example in test/test.py

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



