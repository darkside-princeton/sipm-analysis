#!/usr/bin/env python3
"""Setup Python Modules for Princeton DarkSide Group."""

from setuptools import setup, find_packages
import matplotlib
import os
import shutil


with open('README.md') as f:
  readme = f.read()

with open('requirements.txt') as f:
  required = f.read()

def install_mplstyle():
    stylefile = "darkside.mplstyle"

    mpl_stylelib_dir = os.path.join(matplotlib.get_data_path() ,"stylelib")
    if not os.path.exists(mpl_stylelib_dir):
        os.makedirs(mpl_stylelib_dir)

    print("Installing style into", mpl_stylelib_dir)
    shutil.copy(
        os.path.join(os.path.dirname(__file__), stylefile),
        os.path.join(mpl_stylelib_dir, stylefile))


setup(
  name='sipm',
  version='1.0.0',
  description=('DarkSide-20k Princeton Group Python Analysis Code'),
  long_description=readme,
  author='Ako Jamil',
  author_email='akojamil93@gmail.com',
  url='https://github.com/darkside-princeton/sipm-analysis',
  packages=find_packages(),
  install_requires=required,
  include_package_data=True,
  zip_safe=False,
  package_data= {'styles':['style.mplstyle']},
)

install_mplstyle()
