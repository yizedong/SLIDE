# LSST DECam Image Subtraction Package

This Python package is for performing image subtraction on LSST data with DECam templates.

## Installation

### Prerequisites

Most dependencies of this package has been installed on the Rubin Science Platform. If you miss any packages, you can install them as following:

```bash
pip install --user reproject reproject 
```

### PyZOGY

This package depends on PyZOGY for image subtraction (https://github.com/dguevel/PyZOGY/tree/master):

```bash
# Install PyZOGY from the local directory
cd PyZOGY
pip install --user -e .
```

## Usage

Please see the example.ipynb for how to use this package.



## Citation

If you use this package in your research, please cite: Dong et al. in prep
