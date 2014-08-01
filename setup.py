#!/usr/bin/env python
"""
setup.py file for SWIG kimservice
"""

from __future__ import print_function
from numpy.distutils.core import setup, Extension
from distutils.spawn import find_executable
from subprocess import check_output
import os
import sys

try:
    import configparser
except:
    import ConfigParser as configparser

def run_kim_build_config(opt='--version'):
    buildconf = find_executable("kim-api-build-config")
    if not buildconf:
        print("warning: kim-api-build-config not found on PATH while search for '%s'" % opt)
        return

    return check_output([buildconf, opt])

def getdefault(conf, sect, key, default=None):
    try:
        return conf.get(sect, key)
    except configparser.NoSectionError as e:
        return default
    except configparser.NoOptionError as e:
        return default

def default_kim_dir():
    default = '/usr/local/lib/kim-api/'
    return (
        os.environ.get("KIM_DIR", None) or
        (default if os.path.exists(default) else '') or
        os.path.join(os.environ.get("HOME", None), 'openkim-api')
    )

def default_include():
    libdir = (
        run_kim_build_config('--includes') or
        '-I'+os.path.join(default_kim_dir(), 'include') or ''
    )
    return [ a.replace('-I', '') for a in libdir.split(' ') if a.startswith('-I') ]

def default_libdir():
    libdir = (
        run_kim_build_config('--ldflags') or
        '-L'+os.path.join(default_kim_dir(), '') or ''
    )
    return [ a.replace('-L', '') for a in libdir.split(' ') if a.startswith('-L') ]

def default_libs():
    liblist = run_kim_build_config('--ldlibs') or "-lkim -lgfortran"
    return [ a.strip() for a in liblist.replace('-l', '').split(' ') ]

options = {}
conf = configparser.SafeConfigParser()
setup_cfg = 'setup.cfg' if os.path.exists('setup.cfg') else 'setup.cfg.example'
conf.read(setup_cfg)

options['include'] = getdefault(conf, 'kim-api', 'include') or default_include()
options['libdir'] = getdefault(conf, 'kim-api', 'libdir') or default_libdir()
options['libs'] = getdefault(conf, 'kim-api', 'libraries') or default_libs()

for k,v in options.iteritems():
    if not isinstance(v, list):
        v = [a.strip() for a in v.split(' ')]
        options[k] = v

kimservice_module = Extension('_kimservice',
    sources=['kimservice.i'],
    swig_opts=['-c++'],
    include_dirs=options['include'],
    library_dirs=options['libdir'],
    libraries=options['libs']
)

kimneighborlist_module = Extension('_kimneighborlist',
    sources=['neighborlist.i','neighborlist.c','cvec.c'],
    include_dirs=options['include']+['.'],
    library_dirs=options['libdir'],
    libraries=options['libs']
)

setup(name = 'kimservice',
    version = '0.9.0',
    author      = "Woosong Choi, Matt Bierbaum, Yanjiun Chen",
    description = """KIM python interface""",
    ext_modules = [kimservice_module, kimneighborlist_module],
    py_modules = ['kimdescriptor'],
    data_files=[('kim_virial_template', ['virial/virial.c', 'virial/virial.h', 'virial/virial.i','virial/setup.py'])],
)
