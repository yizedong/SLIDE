# SLIDE: Subtracting LSST Images with DECam Exposures

SLIDE performs image subtraction on LSST data using DECam templates. It is designed to run directly on the Rubin Science Platform and is based on the PyZOGY tutorial by Griffin Hosseinzadeh: https://github.com/griffin-h/image_subtraction

## Installation

### Prerequisites

Most dependencies of this package has been installed on the Rubin Science Platform. If you miss any packages, you can install them as following:

```bash
pip install --user reproject
```

### PyZOGY

This package depends on PyZOGY for image subtraction (https://github.com/dguevel/PyZOGY/tree/master):

```bash
# Install PyZOGY from the local directory
git clone https://github.com/dguevel/PyZOGY.git
cd PyZOGY
pip install --user -e .
```

### Install lsst_decam_subtraction

```bash
# Install lsst_decam_subtraction from the local directory
git clone https://github.com/yizedong/lsst_decam_subtraction.git
cd lsst_decam_subtraction
pip install --user -e .
```

## Documentation

For detailed usage examples, see the `example.ipynb` notebook included in the package.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this package in your research, please cite:

```
Dong et al. (in prep)
```
