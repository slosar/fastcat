# fastcat

Fast and dirty creation of fast and dirty mock galaxy catalogs

The basic module is now fastcat, that just provides a convenient
interface to store galaxy catalogs, together with objectified
understanding of window function and photo-z errors. It is primarily
use as storege inside (the private) DESC 2pt_validation repo for mocks
generate by CoLoRe.

Directories:
* fastcat is the actual python  module.

* hdf5_format_doc is the description of HDF5 format, per version.


# installing

Please run

```
./setup.py install
```

# History

It used to be able to generate them using randomfield package, but
this is now broken.  See attick directory, although not sure how much
is salvageable.




