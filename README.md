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

## TODO

* Shears not yet implemented.



