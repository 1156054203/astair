#-*- coding: utf-8 -*-

import sys
import unittest
from os import path
from setuptools import setup, find_packages


current = path.abspath(path.dirname(__file__))

with open(path.join(current, 'README.md'), 'r') as readme:
    long_description = readme.read()

exec(open('astair/version.py').read())

setup(
    name="asTair",
    version=__version__,
    packages=find_packages(exclude=['tests']),
    install_requires=['Click >=7, < 8', 'pysam >= 0.15.0', 'pyahocorasick >= 1, < 2', 'numpy >= 1, < 2'],
    extras_require={'plot':  ["matplotlib"],},
    test_suite='tests',
    #scripts=['./astair/safe_division.py', './astair/bam_file_parser.py', './astair/simple_fasta_parser.py', './astair/DNA_sequences_operations.py', './astair/context_search.py', './astair/context_search.py', './astair/statistics_summary.py'],
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, < 4',
    author="Gergana V. Velikova and Benjamin Schuster-Boeckler",
    author_email="gergana_velikova@yahoo.com",
    description="A tool for the analysis of bisulfite-free and base-resolution sequencing data generated with TET Assisted Pyridine borane Sequencing (TAPS), or other modified cytosine to thymine conversion methods (mCtoT). It also has some features for bisulfite sequencing data (unmodified cytosine to thymine conversion methods, CtoT).",            long_description=long_description,
    long_description_content_type='text/markdown',
    license="GPLv3",
    entry_points={
        'console_scripts':
        ['astair=scripts.run:cli']
    }, keywords="TAPS taps cytosine caller methylation modification WGBS RRBS bisulfite epigenetics", url="https://bitbucket.org/bsblabludwig/astair/", classifiers=['Programming Language :: Python :: 2.7', 'Programming Language :: Python :: 3.5', 'Programming Language :: Python :: 3.6', 'Programming Language :: Python :: 3.7', 'Intended Audience :: Science/Research', 'Topic :: Scientific/Engineering :: Bio-Informatics']
        )

