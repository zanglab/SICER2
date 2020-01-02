import os
import sys
import glob
from distutils.core import setup
from setuptools import find_packages, Extension
from setuptools.command.build_ext import build_ext as _build_ext

if (float(sys.version[:3])<3):
    sys.stderr.write("ERROR: Python3 required! \n")
    sys.exit(1)

USE_CYTHON = True

EXT = '.pyx' if USE_CYTHON else '.cpp'

extra_cpp_args = ["-O3","-ffast-math", "-stdlib=libc++"]

extensions = [
            Extension('sicer.shared.data_classes',
                depends=glob.glob('sicer/shared/*.h'),
                include_dirs=['.'],
                sources=['sicer/shared/data_classes' + EXT],
                extra_compile_args=extra_cpp_args,
                language='c++'
                ),
            Extension('sicer.shared.chrom_collections',
                include_dirs=['.'],
                sources=['sicer/shared/chrom_collections' + EXT],
                extra_compile_args=extra_cpp_args,
                language='c++'
                ),
            Extension('sicer.bed_reader',
                sources=['sicer/bed_reader' + EXT],
                include_dirs=['.'],
                extra_compile_args=extra_cpp_args,
                language='c++'
                ),
            Extension('sicer.utils_cpp',
                sources=['sicer/utils_cpp' + EXT],
                include_dirs=['.'],
                extra_compile_args=extra_cpp_args,
                language='c++'
                ),
            Extension('sicer.coarsegraining',
                sources=['sicer/coarsegraining' + EXT],
                include_dirs=['.'],
                extra_compile_args=extra_cpp_args
                )
            ]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions, language_level="3")

print(find_packages)

setup(
    name='SICER2',
    version='1.0.1',
    description='SICER2, a redesigned and improved ChIP-seq broad peak calling tool',
    long_description='Redesigned and improved version of the original ChIP-seq broad peak calling tool SICER. Also contains Coarse-graining Approach for Identifying Broad Domains from ChIP-Enriched Regions (RECOGNICER)',
    url='http://zanglab.github.io/SICER2 ',
    author='Jin Yong Yoo, Yiren Wang, Chongzhi Zang*',
    author_email='zang@virginia.edu',
    license='MIT',
    packages=find_packages(),
    package_data={'sicer': ['*.pyx', '*.pxd', '*.h', '*.cpp', '*.c'],
                'sicer.shared': ['*.pyx', '*.pxd', '*.h', '*.cpp', '*.c']
                },
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
    ext_modules=extensions,
    zip_safe=False
)