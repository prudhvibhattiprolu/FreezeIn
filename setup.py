#This program compiles the Freeze-in code and creates a python package
#To manually build: python setup.py --quiet build_ext --inplace clean --all

######################
# Standard Libraries #
######################

import os
import sys
from setuptools import setup, find_packages

##################
# Pybind Library #
##################

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(DIR, "extern", "pybind11"))
from pybind11.setup_helpers import Pybind11Extension, build_ext
del sys.path[-1]

###################################
# Compiling Freeze-in C++ Library #
###################################

FreezeIn = Pybind11Extension(
    #Generate FreezeIn.*.so file in FreezeIn inside the src folder 
    'FreezeIn.FreezeIn',
    sources=['src/FreezeIn/FreezeInPy.cc'],
    language = 'c++',
    define_macros=[('GSTARPATH', '"' + DIR + '"')],
    include_dirs = [
        os.path.join(DIR, 'extern'),
        os.path.join(DIR, 'extern', 'pybind11', 'include')
    ],
    extra_compile_args = [
        #"-H",
        '-std=c++14'
    ]
)

print(DIR);

####################################
# Creating FreezeIn Python Library #
####################################

with open("README.md", "r") as README:
    long_description = README.read()

setup(
    name="FreezeIn",
    version="2.0",
    author="Prudhvi Bhattiprolu",
    author_email="prudhvibhattiprolu@gmail.com",
    description="Computes the portal coupling that reproduces the observed\
            dark matter relic abundance as a function of its mass, for dark\
            matter frozen-in via a light dark photon mediator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    url="https://github.com/prudhvibhattiprolu/FreezeIn",
    packages=find_packages(where="src"),
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
