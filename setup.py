#!/usr/bin/env python3
"""
Setup script for LSST DECam Image Subtraction Package
"""

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return "LSST DECam Image Subtraction Package"


setup(
    name="lsst-decam-subtraction",
    version="0.1.0",
    author="Yize Dong",
    author_email="yize.dong@outlook.com",
        description="Image subtraction package for LSST DECam data",
    install_requires=[
        'numpy>=1.20.0', 'scipy>=1.7.0', 'matplotlib>=3.3.0',
        'astropy>=5.0.0', 'photutils>=1.4.0', 'reproject>=0.8.0', 'sep>=1.1.0',
        'requests>=2.25.0', 'pyvo>=1.3.0', 'astroquery>=0.4.6', 'Pillow>=8.0.0'
    ],
    packages=['lsst_decam_subtraction'],
    entry_points={'console_scripts': ['lsst-decam-subtract=lsst_decam_subtraction.main:main']}
)