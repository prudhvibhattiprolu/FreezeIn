#This program compiles the Freeze-in code and creates a python package
#To manually build: python setup.py --quiet build_ext --inplace clean --all

######################
# Standard Libraries #
######################

import os
import sys
import platform
from setuptools import setup, find_packages

##################
# Pybind Library #
##################

# Set MACOSX_DEPLOYMENT_TARGET only if we are on macOS
if platform.system() == 'Darwin':
    os.environ['MACOSX_DEPLOYMENT_TARGET'] = '10.15'

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(DIR, "extern", "pybind"))
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
        os.path.join(DIR, 'extern', 'pybind', 'include')
    ],
    extra_compile_args = [
#        "-H",
#        '-std=c++14'
    ]
)

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
#    package_data={
#        '': ['gstar/*.tab']
#    },
#    setup_requires=['setuptools_scm'],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
#    install_requires=["numpy>=1.19.1", "scipy>=1.4.1","mpmath>=1.1.0"],
    cmdclass={"build_ext": build_ext},
    ext_modules=[FreezeIn]
)
