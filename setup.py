"""
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.

Copyright (C) 2025 Gran Telescopio Canarias <https://www.gtc.iac.es>
Fabricio Manuel Pérez Toledo <fabricio.perez@gtc.iac.es>
"""

from setuptools import setup, find_packages

setup(
    name="sausero",
    version="0.2.2",
    packages=find_packages(where='.'),
    package_data={
        'SAUSERO    ': [
            'config/configuration.json',
            'BPM/BPM_OSIRIS_PLUS.fits'
        ],
    },
    include_package_data=True,
    install_requires=[
        "astroalign",
	"astrometry_net_client",
	"astropy",
	"astroquery",
	"ccdproc",
	"lacosmic",
	"loguru",
	"matplotlib",
	"numpy",
	"PyYAML",
	"sep"
    ],
    entry_points={
        'console_scripts': [
            'sausero=SAUSERO.OsirisDRP:run',
        ],
    },
    author="Fabricio M. Pérez-Toledo",
    author_email="fabricio.telescope@gmail.com",
    description="This software is designed to reduce Broad Band Imaging observations obtained with OSIRIS+.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/Kennicutt/SAUSERO",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
)
