# LSST DECam Image Subtraction Package

A Python package for performing image subtraction on LSST DECam data, including PSF modeling, reference image downloading, and difference image analysis.

## Installation

### Prerequisites

This package requires several dependencies. Install them using pip:

```bash
pip install astropy numpy matplotlib scipy sep photutils reproject requests pyvo astroquery
```

### PyZOGY

This package depends on PyZOGY for image subtraction:

```bash
# Install PyZOGY from the local directory
cd PyZOGY
pip install --user -e .
```

## Usage

### Command Line Interface

The package provides a command-line interface for easy use:

```bash
python -m lsst_decam_subtraction.subtraction --targimg --ra 58.335 --dec -48.750 --download_DES_temp
```



## Citation

If you use this package in your research, please cite:
