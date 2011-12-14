#!/usr/bin/env python
"""
setup.py file for SWIG kim_virial

Example callback function for virial calculation
"""

from numpy.distutils.core import setup, Extension
import os

kim_dir = os.environ['KIM_DIR']
kim_api_dir = kim_dir + 'KIM_API/'

virial_module = Extension('_virial',
    sources=['virial.i','virial.c'],
    include_dirs=[kim_api_dir, '.'],
    library_dirs=[kim_api_dir],
    libraries=['kim']
    )

setup (name = 'kim_virial',
    version = '0.0.1a',
    author      = "woosong",
    description = """KIM python virial callback""",
    ext_modules = [virial_module],
    )
