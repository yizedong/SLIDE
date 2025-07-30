"""
SLIDE: Subtracting LSST Images with DECam Exposures

This package provides tools for performing image subtraction on LSST images using DECam templates,
including PSF modeling, reference image downloading, and image subtraction.
"""

__version__ = "0.1.0"
__author__ = "Yize Dong"

from .lsst_utils import query_lsst_visits, lsst_visit_to_ccddata, lsst_visit_to_psf
from .image_processing import (
    read_with_datasec,
    get_ccd_bbox,
    make_psf,
    refine_wcs,
    assemble_reference
)
from .reference_download import (
    download_des_reference,
    download_decals_reference,
    gaia3cat,
)
from .subtraction import perform_image_subtraction, lsst_decam_data_load, load_usesr_decam

__all__ = [
    'gaia3cat',
    'query_lsst_visits',
    'read_with_datasec',
    'get_ccd_bbox',
    'make_psf',
    'refine_wcs',
    'assemble_reference',
    'download_des_reference',
    'download_decals_reference',
    'lsst_visit_to_ccddata',
    'lsst_visit_to_psf',
    'perform_image_subtraction',
    'lsst_decam_data_load',
    'load_usesr_decam'
] 