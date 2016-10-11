#!/usr/bin/env python
import distutils
from distutils.core import setup

description = "Fast and dirty cosmological catalogs."

setup(name="fastcat", 
      version="0.9",
      description=description,
      url="https://github.com/slosar/fastcat",
      author="Anze Slosar for LSST DESC",
      author_email="anze@bnl.gov",
      packages=['fastcat'])
