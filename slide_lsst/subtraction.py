"""
Main image subtraction functionality for LSST DECam data.
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sep
import gc
from astropy.nddata import CCDData
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs.utils import skycoord_to_pixel
from astropy.visualization import PercentileInterval, ImageNormalize
from astropy.nddata import Cutout2D
import logging
from PyZOGY.subtract import (
    calculate_difference_image, 
    calculate_difference_image_zero_point,
    normalize_difference_image, 
    save_difference_image_to_file, 
    calculate_difference_psf, 
    save_difference_psf_to_file
)
from PyZOGY.image_class import ImageClass

from .image_processing import (
    read_with_datasec, 
    make_psf, 
    refine_wcs, 
    assemble_reference,
    refine_wcs_astropy,
)
from .reference_download import (
    download_des_reference,
    download_decals_reference,
    gaia3cat,
    des2cat
)
from .lsst_utils import query_lsst_visits, lsst_pixel_to_world, lsst_world_to_pixel, astropy_world_to_pixel, astropy_pixel_to_world, lsst_visit_to_ccddata, lsst_visit_to_psf, lsst_cutout_to_ccddata, safe_cutout2d, lsst_visit_to_psf_median, get_visit_fwhm
from astropy.nddata import CCDData
from astropy import log
import gc 
import psutil
import os
import time
logger = logging.getLogger("slide_lsst")  # or __name__
logger.setLevel(logging.INFO)

des_pixel_scale = 0.262
des_avg_fwhm = {'g':1.11, 
                'r': 0.95, 
                'i': 0.88, 
                'z': 0.83}



def get_memory_usage_gb_cgroup2():
    try:
        with open("/sys/fs/cgroup/memory.current") as f:
            usage_bytes = int(f.read())
        return usage_bytes / 1e9  # Convert to GB
    except FileNotFoundError:
        return None

def wait_for_cgroup_memory(threshold_gb=14.0, check_interval=10):
    """
    Wait until container memory usage (cgroup v2) is below a given threshold.
    
    Parameters
    ----------
    threshold_gb : float
        Max allowed memory usage in GB before pausing.
    check_interval : int
        Seconds to wait before checking again.
    """
    while True:
        mem_gb = get_memory_usage_gb_cgroup2()
        if mem_gb is None:
            print("[Memory Watch] Could not determine memory usage. Proceeding without waiting.")
            break
        if mem_gb < threshold_gb:
            break
        print(f"[Memory Watch] Memory usage is {mem_gb:.2f} GB, waiting to drop below {threshold_gb:.2f} GB...")
        time.sleep(check_interval)

def load_usesr_decam(data, wcs, mask, saturation=65535):
    ccddata = CCDData(data, wcs=wcs, unit='adu', mask=mask)
    ccddata.meta['SATURATE'] = saturation
    return ccddata


def cut_diff_psf(difference_psf):
    real_part = np.real(difference_psf)
    center = np.array(real_part.shape) / 2
    centered_psf = np.roll(real_part, center.astype(int), (0, 1))
    
    #ny, nx = centered_psf.shape
    #center = (nx // 2, ny // 2)
    #cutout = Cutout2D(centered_psf, position=center, size=(25, 25))
    return centered_psf#cutout.data

def save_difference_image_lsst(difference_image, sci_header_primary, sci_header_wcs, normalization, output):
    """
    Save difference image to file.
    
    Normalize and save difference image to file. This also copies over the 
    FITS header of the science image.

    Parameters
    ----------
    difference_image : numpy.ndarray
        Difference image
    science : PyZOGY.ImageClass
        ImageClass instance created from the science image.
    normalization : str
        Normalize to 'reference', 'science', or 'none'.
    output : str
        File to save FITS image to.
    """

    primary_hdu = fits.PrimaryHDU()
    #hdu.header = science.header.copy()
    primary_hdu.header.update(sci_header_primary)
    image_hdu = fits.ImageHDU(data=difference_image, header=sci_header_wcs)

    image_hdu.header['PHOTNORM'] = normalization
    # Combine into an HDUList
    hdul = fits.HDUList([primary_hdu, image_hdu])
    hdul.writeto(output, output_verify='warn', overwrite=True)
    logging.info('Wrote difference image to {}'.format(output))


def perform_image_subtraction(scidata, refdata, sci_psf, ref_psf, show=False, workdir='./', science_filename=None, save_diff=False, sigma_cut=5, max_iterations=3, gain_ratio = np.inf, percent=99, use_pixels=False, size_cut=True, return_output=True, protect_mem=True):
    """
    Perform image subtraction on LSST DECam data.
    
    Parameters
    ----------
    science_filename : str
        Name of the science image file
    template_filename : str, optional
        Name of the template image file (if not downloading)
    workdir : str, optional
        Working directory (default: './')
    ra : float, optional
        Right ascension for reference image query
    dec : float, optional
        Declination for reference image query
    filt : str, optional
        Filter band
    sigma_cut : float, optional
        Sigma cut for difference imaging (default: 5)
    show : bool, optional
        Show results (default: False)
    download_DES_temp : bool, optional
        Download DES template (default: False)
    make_psf_temp : bool, optional
        Make PSF for template (default: False)
        
    Returns
    -------
    str
        Path to output difference image
    """

    if save_diff and science_filename is None:
        science_filename = scidata.meta['filename']

    refdata_aligned = refdata
    
    # Calculate background
    #try:
    #    _data = scidata.data.astype(np.float32)
    #except:
    #    _data = scidata.data.byteswap().newbyteorder()
    #bkg = sep.Background(_data)
    sci_global_bkg = 0 #bkg.globalback
    
    # Perform image subtraction
    logger.info('loading the science image')
    science = ImageClass(scidata, sci_psf.data, scidata.mask, saturation=scidata.meta['SATURATE'])
    science = np.nan_to_num(science, nan=sci_global_bkg) #visit images are already bkg subtracted


    refdata_aligned.mask[np.isnan(refdata_aligned.data)] = True
    logger.info('loading the template image')
    reference = ImageClass(refdata_aligned, ref_psf.data, refdata_aligned.mask, saturation=refdata_aligned.meta['SATURATE'])
    reference = np.nan_to_num(reference, nan=refdata_aligned.meta['ref_global_bkg'])
    combined_mask = reference.mask | scidata.mask
    output_wcs = scidata.wcs  # Save before deleting
    if not return_output:
        del scidata, refdata_aligned
        gc.collect()

    logger.info('calculating the difference image')
    logger.info(f'current memory usage: {get_memory_usage_gb_cgroup2()}')
    if protect_mem:
        wait_for_cgroup_memory(threshold_gb=15, check_interval=15)
    difference = calculate_difference_image(science, reference, show=show, max_iterations=max_iterations, sigma_cut=sigma_cut, gain_ratio=gain_ratio, percent=percent, use_pixels=use_pixels, size_cut=size_cut)
    difference_zero_point = calculate_difference_image_zero_point(science, reference)
    normalized_difference = normalize_difference_image(difference, difference_zero_point, science, reference, 'i')
    del science, reference, difference
    gc.collect()
    logger.info(f"Image subtraction completed successfully!")
    if save_diff:
        output_filename = os.path.join(workdir, science_filename.replace('.fits', '.diff.fits'))
        primary_hdu = fits.PrimaryHDU(normalized_difference, header=output_wcs.to_header())
        mask_data = np.array(combined_mask, dtype=np.uint8)
        mask_hdu = fits.ImageHDU(mask_data, name='MASK')
        psf_data = sci_psf.data #cutout_difference_psf#
        psf_hdu = fits.ImageHDU(psf_data, name='PSF')
        hdul = fits.HDUList([primary_hdu, mask_hdu, psf_hdu])
        hdul.writeto(output_filename, overwrite=True)
    if return_output:
        normalized_difference = CCDData(normalized_difference, wcs=output_wcs, unit='adu')
        #normalized_difference = np.where(combined_mask, np.nan, normalized_difference)
        normalized_difference.mask = combined_mask
        return normalized_difference, sci_psf.data
    else:
        return None

def lsst_decam_data_load(visit_image, ra=None, dec=None, science_filename = None, template_filename=None, user_decam_data=None, workdir='./', show=False, download_DES_temp=False, download_DECaLS_temp=False, cutout=False, cutout_size=1000, get_median_sci_psf=True, make_sci_psf=False, reference_catalog='gaia', reference_mag1=17, reference_mag2=21, save_intermediate=False, save_original_temp=False, fit_distortion=None, refine_wcs_sci=False,
                        refine_wcs_ref=False, mask_type = []):
    """
    Perform image subtraction on LSST DECam data.
    """
    if len(mask_type) == 0:
        mask_type = ["BAD", "SAT", "INTRP", "CR", "EDGE", "SUSPECT", "NO_DATA", "SENSOR_EDGE", "CLIPPED", "CROSSTALK", "UNMASKEDNAN", "STREAK"]
    image_filter = visit_image.filter.bandLabel

    if science_filename is None:
        visit_info = visit_image.getInfo().getVisitInfo()
        visit_id = visit_info.getId()
        detector_id = visit_image.getDetector().getId()
        science_filename = f'{visit_id}_{detector_id}_{image_filter}.fits'
    
    # Read the science image
    if cutout:
        scidata = safe_cutout2d(visit_image, ra, dec, cutout_size=cutout_size, mask_type = mask_type)
    else:
        scidata = lsst_visit_to_ccddata(visit_image, mask_type = mask_type)
        ra, dec = lsst_pixel_to_world(2036.0, 2000.0, visit_image)
    logger.info(f'science image saturation level: {scidata.meta['SATURATE']}')

    nx, ny = scidata.shape
    half_ra, half_dec = astropy_pixel_to_world(nx//2, ny//2, scidata.wcs)

    # Get Gaia/DES catalog for PSF modeling
    if reference_catalog == 'gaia':
        catalog = gaia3cat(ra=np.round(half_ra, 3), dec=np.round(half_dec, 3), ccddata=scidata, band=image_filter, radius_arcmin=11, mag1=reference_mag1, mag2=reference_mag2)
        des_fwhm = des_avg_fwhm[image_filter]/des_pixel_scale
    elif reference_catalog == 'des': 
        catalog = des2cat(ra=np.round(half_ra, 3), dec=np.round(half_dec, 3), ccddata=scidata, band=image_filter, radius_arcmin=11, mag1=reference_mag1, mag2=reference_mag2)
        des_fwhm = np.median(catalog[f'{image_filter}FWHM'])
    logger.info(f'catalog size: {len(catalog)}')
    catalog['raMean'], catalog['decMean'] = catalog['ra'], catalog['dec']

    if refine_wcs_sci:
    # Refine WCS for science image
        logger.info(f"Using {len(catalog)} stars for WCS refinement")
        new_wcs, catalog = refine_wcs_astropy(scidata.data, scidata.wcs, catalog, fwhm=get_visit_fwhm(visit_image, ra, dec), fit_distortion=fit_distortion)
        scidata.wcs = new_wcs

    if save_intermediate:
        _science_filename = os.path.join(workdir, science_filename)
        scidata.write(_science_filename, overwrite=True)
    
    # Read science image PSF
    psf_flag = 0
    if make_sci_psf:
        if len(catalog) < 5:
            logger.info('Too few stars, getting PSF from the visit_image')
            sci_psf = lsst_visit_to_psf_median(visit_image, ra, dec, cutout_size=cutout_size)
            psf_flag=1
        else:
            sci_psf, _ = make_psf(scidata, catalog, show=show)
    elif get_median_sci_psf:
        sci_psf = lsst_visit_to_psf_median(visit_image, ra, dec, cutout_size=cutout_size)
        psf_flag = 1
    else:
        sci_psf = lsst_visit_to_psf(visit_image, ra, dec)
        psf_flag = 1
    if show and psf_flag==1:
        fig, ax1 = plt.subplots(1, 1, figsize=(5, 5))
        norm = ImageNormalize(sci_psf, PercentileInterval(98))
        ax1.imshow(sci_psf, norm=norm)
        ax1.set_title('Science PSF')

    # Get reference image
    if download_DES_temp:
        if ra is None or dec is None:
            raise ValueError("RA and DEC must be provided when downloading DES template")
        _refdata = download_des_reference(ra=half_ra, dec=half_dec, fov=max(nx, ny)*0.2/3600*1.5, filt=image_filter)
        logger.info(f'DES template downloaded at RA:{ra} DEC:{dec} in the {image_filter} band')
        if _refdata == None:
            logger.info(f'no DES images found; attempting to download decals templates at RA:{ra} DEC:{dec} in the {image_filter} band')
            _refdata = download_decals_reference(ra=half_ra, dec=half_dec, fov=max(nx, ny)*0.2/3600*1.5, filt=image_filter)
            logger.info(f'decals template downloaded at RA:{ra} DEC:{dec} in the {image_filter} band')
        if save_original_temp:
            _science_filename = os.path.join(workdir, science_filename.replace('.fits', '.oritemp.fits'))
            _refdata.write(_science_filename, overwrite=True)
    elif download_DECaLS_temp:
        if ra is None or dec is None:
            raise ValueError("RA and DEC must be provided when downloading decals template")
        _refdata = download_decals_reference(ra=half_ra, dec=half_dec, fov=max(nx, ny)*0.2/3600*1.5, filt=image_filter)
        if save_original_temp:
            _science_filename = os.path.join(workdir, science_filename.replace('.fits', '.oritemp.fits'))
            _refdata.write(_science_filename, overwrite=True)
        logger.info(f'decals template downloaded at RA:{ra} DEC:{dec} in the {image_filter} band')
    else:
        if user_decam_data is None:
            raise ValueError("A DeCam image needs to be provided")
        _refdata = user_decam_data
        logger.info(f'using user-provided decam template')

    if refine_wcs_ref:
        # Refine WCS for reference image
        new_wcs, catalog = refine_wcs_astropy(_refdata.data, _refdata.wcs, catalog, fwhm=des_fwhm, fit_distortion=fit_distortion)
        _refdata.wcs = new_wcs

    # Calculate background
    try:
        _data = _refdata.data.astype(np.float32)
    except:
        _data = _refdata.data.byteswap().newbyteorder()
    bkg = sep.Background(_data)
    ref_global_bkg = bkg.globalback

    # Assemble reference image
    logger.info('align the template with the science image')
    refdata_aligned = assemble_reference([_refdata], scidata.wcs, scidata.shape, ref_global_bkg=ref_global_bkg)
    refdata_aligned.meta['SATURATE'] = _refdata.meta['SATURATE']
    if save_intermediate:
        _template_filename = os.path.join(workdir, science_filename.replace('.fits', '.temp.fits'))
        refdata_aligned.write(_template_filename, overwrite=True)
        
    # Get reference PSF
    logger.info(f'Making PSF for the template image')
    ref_psf, _ = make_psf(refdata_aligned, catalog, show=show)

    scidata.meta['filename'] = science_filename
    refdata_aligned.meta['ref_global_bkg'] = ref_global_bkg
    return scidata, refdata_aligned, sci_psf, ref_psf

    