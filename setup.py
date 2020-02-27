import os
import sys
import warnings
import subprocess
from distutils.core import setup
from setuptools import find_packages, Extension
from setuptools.command.build_ext import build_ext as _build_ext

if (float(sys.version[:3])<3):
    sys.stderr.write('ERROR: Python3 required! \n')
    sys.exit(1)

if (float(sys.version[:3])<3.5):
    def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
        return '%s:%s: %s:%s\n' % (filename, lineno, category.__name__, message)

    warnings.formatwarning = warning_on_one_line
    warnings.warn("Recommended to use Python 3.5 or above to run SICER2.")

USE_CYTHON = False

EXT = '.pyx' if USE_CYTHON else '.cpp'

extra_cpp_args = ['-O3','-ffast-math', '-w', '-std=c++11']

# Check if C++ compiler is GCC or Clang
if "clang" in subprocess.check_output(['gcc', '--version']).decode('utf-8').lower():
    extra_cpp_args.append('-stdlib=libc++')

extension_names = [
    'sicer.shared.data_classes', 'sicer.shared.containers', 'sicer.shared.utils',
    'sicer.file_writers', 'sicer.bed_reader', 'sicer.generate_windows', 
    'sicer.find_islands', 'sicer.associate_tags_with_control', 
    'sicer.filter_islands_by_fdr', 'sicer.recover_significant_reads', 
    'sicer.coarsegraining', 'sicer.find_union_islands', 'sicer.compare_two_libraries'
    ]

def generate_extensions(name_list):
    extensions = []
    for name in name_list:
        extension = Extension(
            name, 
            include_dirs=['.', 'sicer/shared'],
            sources=[name.replace('.', '/') + EXT],
            extra_compile_args=extra_cpp_args,
            language='c++'
        )
        extensions.append(extension)

    return extensions


extensions = generate_extensions(extension_names)

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions, language_level='3')

data_ext = ['*.pyx', '*.pxd', '*.h', '*.c', '*.hpp', '*.cpp']

setup(
    name='SICER2',
    version='1.1.0',
    description='SICER2, a redesigned and improved ChIP-seq broad peak calling tool',
    long_description='Redesigned and improved version of the original ChIP-seq broad peak calling tool SICER. Also contains Coarse-graining Approach for Identifying Broad Domains from ChIP-Enriched Regions (RECOGNICER)',
    url='http://zanglab.github.io/SICER2 ',
    author='Jin Yong Yoo, Yiren Wang, Chongzhi Zang*',
    author_email='zang@virginia.edu',
    license='MIT',
    packages=find_packages(),
    package_data={'sicer': data_ext + ['genomedata/*.json'], 'sicer.shared': data_ext},
    scripts=['bin/sicer','bin/sicer_df', 'bin/recognicer', 'bin/recognicer_df'],
    setup_requires=['numpy','scipy>=1.0.0'],
    install_requires=['numpy','scipy>=1.0.0'],
    keywords = ['ChIP-Seq','SICER'],
    classifiers=['Programming Language :: Python :: 3',
        'Environment :: Other Environment',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
        'Topic :: Scientific/Engineering'],
    ext_modules=extensions,
    zip_safe=False
)