#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import re

try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup

setup(
    name="mkc3kgrid",
    version='0.1',
    author="Phill Cargile",
    author_email="pcargile@cfa.harvard.edu",
    packages=["mkc3kgrid"],
    url="",
    #license="LICENSE",
    description="Generate grid of spectra/phot models directly from Payne ANN",
    long_description=open("README.md").read() + "\n\n",
    #install_requires=["numpy", "scipy"],
)
