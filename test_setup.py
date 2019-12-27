from distutils.core import setup
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

extensions = [
        Extension('sicer.src.read_bed',
                  sources=['sicer/src/read_bed.pyx'],
                  extra_compile_args=['-O3', '-stdlib=libc++'],
                  language='c++'),
        Extension('sicer.src.shared.read',
                  sources=['sicer/src/shared/read.pyx'],
                  extra_compile_args=['-O3', '-stdlib=libc++'],
                  language='c++')
        ]

setup(ext_modules=cythonize(extensions),
    package_data = {'sicer': ['sicer/src/*.pxd']}
    )
