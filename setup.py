from setuptools import setup, find_packages
from os import path

current = path.abspath(path.dirname(__file__))

with open(path.join(current, 'README.md'), encoding='utf-8') as readme:
long_description = readme.read()

setup(
    name="asTair",
    version="3.0",
    packages=find_packages(),
    scripts=['/astair/astair_aligner_v3.py', '/astair/astair_mbias_v3.py', '/astair/astair_mod_caller_v3.py', '/astair/astair_phred_values_v3.py', '/astair/astair_taps_modification_simulation_v3.py'],
    install_requires=['click', 'pysam >= 0.15.0', 'pyahocorasick', 'numpy'],
    extras_require={'Plot':  ["matplotlib"],},
    python_requires='>=2.7', '>=3.5',
    author="Gergana V. Velikova and Benjamin Schuster-Boeckler",
    author_email="gergana_velikova@yahoo.com",
    description="A little tool developed for the analysis of bisulfite-free and base-resolution sequencing data generated with TET Assisted Pyridine borane Sequencing (TAPS) or other modified cytosine to thymine conversion methods (mCtoT), but it also has some tricks for bisulfite sequencing data (unmodified cytosine to thymine conversion methods, CtoT).",
    long_description=long_description,
    long_description_content_type='text/markdown'
    license="GPLv3",
    keywords="TAPS taps cytosine caller methylation modification WGBS RRBS bisulfite epigenetics",
    url="https://bitbucket.org/bsblabludwig/astair/"
    classifiers=['Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'])

