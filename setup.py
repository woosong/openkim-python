#!/usr/bin/env python
"""
setup.py file for SWIG kimservice
"""

from numpy.distutils.core import setup, Extension
import os

kim_dir = os.environ['KIM_DIR']
kim_api_dir = kim_dir + 'KIM_API/'

kimservice_module = Extension('_kimservice',
    sources=['kimservice.i'],
    swig_opts=['-c++'],
    include_dirs=[kim_api_dir],
    library_dirs=[kim_api_dir],
    libraries=['kim']
    )

neighborlist_module = Extension('_neighborlist',
    sources=['neighborlist.i','neighborlist.c'],
    include_dirs=[kim_api_dir, '.'],
    library_dirs=[kim_api_dir],
    libraries=['kim']
    )

setup (name = 'kimservice',
    version = '0.0.1a',
    author      = "woosong",
    description = """KIM python interface""",
    ext_modules = [kimservice_module, neighborlist_module],
    py_modules = ['kimdescriptor'],
    data_files=[('kim_virial_template', ['virial/virial.c', 'virial/virial.h', 'virial/virial.i','virial/setup.py'])],
    )
