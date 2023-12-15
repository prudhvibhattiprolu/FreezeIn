#To build:
# >>> python setup.py --quiet build_ext --inplace clean --all

import os
import sys
from glob import glob
import setuptools
from setuptools import setup

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(DIR, "extern", "pybind11"))

from pybind11.setup_helpers import Pybind11Extension, build_ext

del sys.path[-1]

with open("README.md", "r") as README:
    long_description = README.read()

FreezeIn = Pybind11Extension(
    'FreezeIn',
    sources=['src/FreezeInPy.cc'],
    language = 'c++',
    include_dirs = [
        os.path.join(DIR, 'extern'),
        os.path.join(DIR, 'extern', 'pybind11', 'include')
    ],
    extra_compile_args = [
        #"-H",
        '-std=c++14'
    ]
)

setup(
    name="FreezeIn",
    version="1.0",
    author="Prudhvi Bhattiprolu",
    author_email="prudhvibhattiprolu@gmail.com",
    description="Computes the direct detection cross section or equivalently\
            the portal couplings that reproduce the observed relic abundance,\
            as a function of its mass, for dark matter frozen-in via a light\
            dark photon mediator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    url="https://github.com/prudhvibhattiprolu/FreezeIn",
    packages=setuptools.find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    cmdclass={"build_ext": build_ext},
    ext_modules=[FreezeIn]
)
