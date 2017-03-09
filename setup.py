#! /usr/bin/env python
#
# Copyright (C) 2016 Ian Driver <ian.driver@ucsf.edu>
import os


DESCRIPTION = "scicast: Single Cell Iterative Clustering and Statistical Testing. A package for interrogating single cell sequencing data."
LONG_DESCRIPTION = """\
scicast is a python utility that automates many of the repetitive steps of analyzing single cell sequencing data.

-k-means clustering to identify clusters

-Clustering and subclustering of data to identify 'stable' sets of cells.

-Statistical testing to identify top genes that indentify stable cluster.

-Correlation search and analysis to identify gene networks driving cluster identity.

-Outputs both plots for visualization (PCA and heatmap) cell and gene lists that can be used to refine analysis.
"""

DISTNAME = 'scicast'
MAINTAINER = 'Ian Driver'
MAINTAINER_EMAIL = 'ian.driver@ucsf.edu'
URL = 'https://github.com/iandriver/scicast'
LICENSE = 'MIT'
DOWNLOAD_URL = 'https://github.com/iandriver/scicast'
VERSION = '0.8.27'

try:
    from setuptools import setup
    _has_setuptools = True
except ImportError:
    from distutils.core import setup

def check_dependencies():
    install_requires = []

    # Just make sure dependencies exist, I haven't rigorously
    # tested what the minimal versions that will work are
    # (help on that would be awesome)
    try:
        import numpy
    except ImportError:
        install_requires.append('numpy')
    try:
        import scipy
    except ImportError:
        install_requires.append('scipy')
    try:
        import sklearn
    except ImportError:
        install_requires.append('scikit-learn')
    try:
        import matplotlib
    except ImportError:
        install_requires.append('matplotlib')
    try:
        import pandas
    except ImportError:
        install_requires.append('pandas>=0.19.0')
    try:
        import seaborn
    except ImportError:
        install_requires.append('seaborn>=0.7.1')
    try:
        import rpy2
    except ImportError:
        install_requires.append('rpy2')
    try:
        import fastcluster
    except ImportError:
        install_requires.append('fastcluster')

    return install_requires

if __name__ == "__main__":

    install_requires = check_dependencies()

    setup(name=DISTNAME,
        author=MAINTAINER,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        license=LICENSE,
        url=URL,
        version=VERSION,
        py_modules=['scicast.cluster', 'scicast.matrix_filter'],
        entry_points={
          'console_scripts': ['scicast = scicast.cluster:main']
        },
        download_url=DOWNLOAD_URL,
        install_requires=install_requires,
        packages=['scicast'],
        keywords='single-cell single cell RNA-seq sequencing clustering PCA k-means',
        classifiers=[
                     'Intended Audience :: Science/Research',
                     'Programming Language :: Python :: 2.7',
                     'Programming Language :: Python :: 3.5',
                     'Programming Language :: Python :: 3.6',
                     'License :: OSI Approved :: MIT License',
                     'Topic :: Scientific/Engineering :: Visualization',
                     'Topic :: Scientific/Engineering :: Bio-Informatics',
                     'Operating System :: POSIX',
                     'Operating System :: Unix',
                     'Operating System :: MacOS']
          )
