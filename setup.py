#!/usr/bin/env python3
"""
Setup script for SLIDE: Subtracting LSST Images with DECam Exposures
"""

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return "SLIDE: Subtracting LSST Images with DECam Exposures"

setup(
    name="slide_lsst",
    version="0.1.0",
    author="Yize Dong",
    author_email="yize.dong@outlook.com",
    description="Image subtraction package for LSST DECam data",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/yizedong/SLIDE",
    project_urls={
        "Bug Tracker": "https://github.com/yizedong/SLIDE/issues",
        "Documentation": "https://github.com/yizedong/SLIDE#readme",
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords="astronomy, lsst, decam, image-subtraction, transient-detection",
    packages=find_packages(),
    python_requires='>=3.10',
    install_requires=[
        'numpy>=1.20.0',
        'scipy>=1.7.0', 
        'matplotlib>=3.3.0',
        'astropy>=5.0.0',
        'photutils>=1.4.0',
        'reproject>=0.8.0',
        'sep>=1.1.0',
        'requests>=2.25.0',
        'pyvo>=1.3.0',
        'astroquery>=0.4.6',
        'Pillow>=8.0.0'
    ],
    extras_require={
        'dev': [
            'pytest>=6.0',
            'pytest-cov>=2.0',
            'black>=21.0',
            'flake8>=3.8',
        ],
        'docs': [
            'sphinx>=4.0',
            'sphinx-rtd-theme>=1.0',
        ],
    },
    include_package_data=True,
    zip_safe=False,
)