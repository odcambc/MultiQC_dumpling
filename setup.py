#!/usr/bin/env python
"""
MultiQC_dumpling is a plugin for MultiQC, providing additional tools which are
specific to Dumpling DMS analysis pipeline.

For more information about Dumpling, see https://github.com/odcambc/dumpling
For more information about MultiQC, see http://multiqc.info
"""

from setuptools import setup, find_packages

version = '0.1'

setup(
    name = 'multiqc_dumpling',
    version = version,
    author = 'Christian Macdonald',
    author_email = 'christian.macdonald@ucsf.edu',
    description = "MultiQC plugin for Dumpling DMS pipeline",
    long_description = __doc__,
    keywords = 'bioinformatics',
    url = 'https://github.com/odcambc/MultiQC_dumpling',
    download_url = 'https://github.com/odcambc/MultiQC_dumpling/releases',
    license = 'MIT',
    packages = find_packages(),
    include_package_data = True,
    zip_safe=False,
    install_requires = [
        'multiqc>=1.0',
    ],
    entry_points = {
        'multiqc.modules.v1': [
            'dumpling = multiqc_dumpling.modules.dumpling:MultiqcModule',
        ],
        'multiqc.hooks.v1': [
            'config_loaded = multiqc_dumpling.dumpling_intialization:dumpling_plugin_execution_start'
        ],
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)