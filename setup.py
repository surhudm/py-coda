#!/usr/bin/env python

from distutils.core import setup

setup(name='py-coda',
      version='0.99',
      description='python bindings for MCMC diagnostics provided by the package coda in R',
      author='Surhud More',
      author_email='surhudkicp@gmail.com',
      url='http://www.github.com/surhudm/py-coda',
      license     = ['GPL'],
      py_modules=['py_coda'],
      package_dir = { '':'src'},
     )
