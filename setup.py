import os
import sys
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as _build_ext


if (float(sys.version[:3])<3):
    sys.stderr.write("ERROR: Python3 required! \n")
    sys.exit(1)

extra_c_args = ["-w","-O3","-ffast-math"]

ext_modules = [Extension("sicer.src.coarsegraining",["sicer/src/coarsegraining.c"],extra_compile_args=extra_c_args)]

setup(
    name='SICER2',
    version='1.0.2',
    description = 'SICER2, a redesigned and improved ChIP-seq broad peak calling tool',
    long_description='Redesigned and improved version of the original ChIP-seq broad peak calling tool SICER. Also contains Coarse-graining Approach for Identifying Broad Domains from ChIP-Enriched Regions (RECOGNICER)',
    url = 'http://zanglab.github.io/SICER2 ',
    author = 'Jin Yong Yoo, Yiren Wang, Chongzhi Zang*',
    author_email = 'zang@virginia.edu',
    license = 'MIT',
    packages=find_packages(),
    scripts=['bin/sicer','bin/sicer_df', 'bin/recognicer', 'bin/recognicer_df'],
    setup_requires=['numpy','scipy>=1.0.0'],
    install_requires=['numpy','scipy>=1.0.0'],
    keywords = ['ChIP-Seq','SICER'],
    classifiers=["Programming Language :: Python :: 3",
        "Environment :: Other Environment",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"],
    ext_modules=ext_modules
)
