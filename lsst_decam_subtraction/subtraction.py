"""
Main image subtraction functionality for LSST DECam data.
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sep
import gc
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs.utils import skycoord_to_pixel
from astropy.visualization import PercentileInterval, ImageNormalize
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
    refine_wcs_astropy
)
from .reference_download import (
    download_decam_reference,
    gaia3cat
)
from .lsst_utils import query_lsst_visits, lsst_pixel_to_world, lsst_world_to_pixel, astropy_world_to_pixel, astropy_pixel_to_world, lsst_visit_to_ccddata, lsst_visit_to_psf, lsst_cutout_to_ccddata, safe_cutout2d
from astropy.nddata import CCDData
from astropy import log
logger = logging.getLogger("lsst_decam_subtraction")  # or __name__
logger.setLevel(logging.INFO)




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


def perform_image_subtraction(scidata, refdata, sci_psf, ref_psf, ref_global_bkg, sci_header=None, save_intermediate=False, show=False, workdir='./', science_filename='test.fits', sigma_cut=5, max_iterations=3, gain_ratio = np.inf, percent=99, use_pixels=False, size_cut=True):
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

    if save_intermediate and science_filename is None:
        raise Warning("Science filename is required to save intermediate files")
    """
    # Calculate background
    try:
        _data = refdata.data.astype(np.float32)
    except:
        _data = refdata.data.byteswap().newbyteorder()
    bkg = sep.Background(_data)
    ref_global_bkg = bkg.globalback

    # Assemble reference image
    logger.info('align the template with the science image')
    refdata_aligned = assemble_reference([refdata], scidata.wcs, scidata.shape, ref_global_bkg=ref_global_bkg)
    
    if save_intermediate and science_filename is not None:
        _template_filename = os.path.join(workdir, science_filename.replace('.fits', '.temp.fits'))
        refdata_aligned.write(_template_filename, overwrite=True)
    """
    refdata_aligned = refdata
    # Perform image subtraction
    logger.info('loading the science image')
    science = ImageClass(scidata, sci_psf.data, scidata.mask, saturation=65535)
    #science.background_counts = np.zeros_like(science.background_counts)
    #print(science.background_counts, science.background_std)

    refdata_aligned.mask[np.isnan(refdata_aligned.data)] = True
    logger.info('loading the template image')
    reference = ImageClass(refdata_aligned, ref_psf.data, refdata_aligned.mask, saturation=65535)
    #reference.background_counts = np.zeros_like(reference.background_counts)
    #print(reference.background_counts, reference.background_std)
    #reference.data = np.nan_to_num(reference.data, nan=ref_global_bkg)
    reference = np.nan_to_num(reference, nan=ref_global_bkg)

    logger.info('calculating the difference image')
    difference = calculate_difference_image(science, reference, show=show, max_iterations=max_iterations, sigma_cut=sigma_cut, gain_ratio=gain_ratio, percent=percent, use_pixels=use_pixels, size_cut=size_cut)
    difference_zero_point = calculate_difference_image_zero_point(science, reference)
    normalized_difference = normalize_difference_image(difference, difference_zero_point, science, reference, 'i')
    if save_intermediate and science_filename is not None:
        output_filename = os.path.join(workdir, science_filename.replace('.fits', '.diff.fits'))
        #save_difference_image_lsst(normalized_difference, sci_header[0], sci_header[1], 'i', output_filename)
        hdu = fits.PrimaryHDU(normalized_difference, header=scidata.wcs.to_header())
        hdu.writeto(output_filename, overwrite=True)
    normalized_difference = CCDData(normalized_difference, wcs=scidata.wcs, unit='adu')
    combined_mask = reference.mask | scidata.mask
    normalized_difference = np.where(combined_mask, np.nan, normalized_difference)
    
    logger.info(f"Image subtraction completed successfully!")

    return refdata_aligned, normalized_difference, sci_psf.data

def lsst_decam_data_load(visit_image, ra=None, dec=None, science_filename = 'test.fits', template_filename=None, workdir='./', show=False, download_DES_temp=False, cutout=False,
                        cutout_size=1000, save_intermediate=False, save_original_temp=False,  
                         fit_distortion=None):
    """
    Perform image subtraction on LSST DECam data.
    """
    
    image_filter = visit_image.filter.bandLabel
    # Read the science image
    if cutout:
        _scidata = safe_cutout2d(visit_image, ra, dec, cutout_size=cutout_size)
        #mask = safe_cutout2d(visit_image.mask.array, ra, dec, cutout_size=cutout_size)
        #_scidata = lsst_cutout_to_ccddata(cutout)
    else:
        _scidata = lsst_visit_to_ccddata(visit_image)
    if save_intermediate and science_filename is not None:
            _science_filename = os.path.join(workdir, science_filename.replace('.fits', '.orisci.fits'))
            _scidata.write(_science_filename, overwrite=True)
    nx, ny = _scidata.shape
    half_ra, half_dec = astropy_pixel_to_world(nx//2, ny//2, _scidata.wcs)

    # Get Gaia catalog for PSF modeling
    catalog = gaia3cat(ra=np.round(half_ra, 3), dec=np.round(half_dec, 3), radius_arcmin=11)
    catalog['raMean'], catalog['decMean'] = catalog['ra'], catalog['dec']
    logger.info(f'Gaia catalog size: {len(catalog)}')
    
    # Filter stars within image bounds
    x, y = astropy_world_to_pixel(catalog['ra'], catalog['dec'], _scidata.wcs)
    margin = 25 // 2
    
    mask = (
        (x > margin) & (x < nx - margin) &
        (y > margin) & (y < ny - margin)
    )
    catalog = catalog[mask]

    # Refine WCS for science image
    logger.info(f"Using {len(catalog)} stars for WCS refinement")
    scidata = _scidata.copy()
    new_wcs = refine_wcs_astropy(scidata.data, scidata.wcs, catalog, fit_distortion=fit_distortion)
    scidata.wcs = new_wcs

    if save_intermediate and science_filename is not None:
        _science_filename = os.path.join(workdir, science_filename)
        scidata.write(_science_filename, overwrite=True)

    # Read science image PSF
    sci_psf = lsst_visit_to_psf(visit_image, ra, dec)

    # Get reference image
    if download_DES_temp:
        if ra is None or dec is None:
            raise ValueError("RA and DEC must be provided when downloading DES template")
        _refdata = download_decam_reference(ra=half_ra, dec=half_dec, fov=0.3, filt=image_filter)
        if save_original_temp and science_filename is not None:
            _science_filename = os.path.join(workdir, science_filename.replace('.fits', '.oritemp.fits'))
            _refdata.write(_science_filename, overwrite=True)
        logger.info(f'DES template downloaded for RA:{ra} DEC:{dec} for {image_filter} band')
    else:
        if template_filename is None:
            raise ValueError("Template filename must be provided when not downloading")
        filename = os.path.join(workdir, template_filename)
        _refdata = read_with_datasec(filename)
    
    # Refine WCS for reference image
    new_wcs = refine_wcs_astropy(_refdata.data, _refdata.wcs, catalog, fit_distortion=fit_distortion)
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
    
    if save_intermediate and science_filename is not None:
        _template_filename = os.path.join(workdir, science_filename.replace('.fits', '.temp.fits'))
        refdata_aligned.write(_template_filename, overwrite=True)
    
    # Get reference PSF
    logger.info(f'Making PSF for the template image')
    ref_psf, _ = make_psf(refdata_aligned, catalog, show=show)

    
    return scidata, refdata_aligned, sci_psf, ref_psf, ref_global_bkg



def local_data_load(ra=None, dec=None, science_filename=None, template_filename=None, workdir='./', filt=None, sigma_cut=5, show=False, download_DES_temp=False, 
                     make_psf_temp=False):
    """
    Perform image subtraction on LSST DECam data.
    """
    # Read the science image
    filename = os.path.join(workdir, science_filename)
    _scidata = read_with_datasec(filename)

    nx, ny = _scidata.shape
    half_ra, half_dec = _scidata.wcs.all_pix2world(nx//2, ny//2, 0)
    
    print(f"Science image center: RA={np.round(half_ra, 3)}, DEC={np.round(half_dec, 3)}")
    
    # Get Gaia catalog for PSF modeling
    catalog = gaia3cat(ra=np.round(half_ra, 3), dec=np.round(half_dec, 3), radius_arcmin=10)
    print(f'Gaia catalog size: {len(catalog)}')
    catalog['raMean'], catalog['decMean'] = catalog['ra'], catalog['dec']
    
    # Filter stars within image bounds
    coords = SkyCoord(catalog['raMean'], catalog['decMean'], unit='deg')
    x, y = skycoord_to_pixel(coords, _scidata.wcs, origin=0)
    margin = 25 // 2
    
    mask = (
        (x > margin) & (x < nx - margin) &
        (y > margin) & (y < ny - margin)
    )
    catalog = catalog[mask]

    # Create PSF for science image
    _, sci_stars = make_psf(_scidata, catalog, show=show, boxsize=25)
    scidata = _scidata.copy()
    refine_wcs(scidata.wcs, sci_stars, catalog)
    print(f"Using {len(catalog)} stars for WCS refinement")

    # Read science image PSF
    filename = os.path.join(workdir, science_filename.replace('.fits', '.psf.fits'))
    sci_psf = read_with_datasec(filename)

    # Get reference image
    if download_DES_temp:
        if ra is None or dec is None:
            raise ValueError("RA and DEC must be provided when downloading DES template")
        _refdatas = download_decam_reference(ra=ra, dec=dec, fov=0.5, filt=filt)
        print(f'DES template downloaded for RA:{ra} DEC:{dec}')
    else:
        if template_filename is None:
            raise ValueError("Template filename must be provided when not downloading")
        filename = os.path.join(workdir, template_filename)
        _refdatas = read_with_datasec(filename)
    
    # Get catalog for reference image
    if ra is not None and dec is not None:
        half_ra, half_dec = ra, dec
    else:
        nx, ny = _refdatas.shape
        half_ra, half_dec = _refdatas.wcs.all_pix2world(nx//2, ny//2, 0)
    
    catalog = gaia3cat(ra=half_ra, dec=half_dec, radius_arcmin=10)
    catalog['raMean'], catalog['decMean'] = catalog['ra'], catalog['dec']
    
    # Create PSF for reference image
    _, ref_stars = make_psf(_refdatas, catalog, show=show)
    refdatas = _refdatas.copy()
    
    
    
    refine_wcs(refdatas.wcs, ref_stars, catalog)

    # Get reference PSF
    if make_psf_temp:
        print(f'Making PSF for {template_filename}')
        ref_psf, _ = make_psf(refdatas, catalog, show=False)
    else:
        filename = os.path.join(workdir, template_filename.replace('.fits', '.psf.fits'))
        print(f'Reading PSF for {template_filename}')
        ref_psf = read_with_datasec(filename)
    
    return scidata, refdata, sci_psf, ref_psf





def main():
    """Main command-line interface."""
    description = "LSST DECam Image Subtraction Tool"
    usage = "%(prog)s [options]"
    parser = argparse.ArgumentParser(usage=usage, description=description)
    
    # Required arguments
    parser.add_argument("--targimg", dest="targimg", required=True,
                       help='Name of the science image file')
    
    # Optional arguments
    parser.add_argument("--tempimg", dest="tempimg", default=None,
                       help='Name of the template image file (if not downloading)')
    parser.add_argument("--imgdir", dest="imgdir", default='./',
                       help='Path to the images directory (default: ./)')
    parser.add_argument("--ra", dest="ra", type=float, default=None,
                       help='RA of the object (used for querying reference image)')
    parser.add_argument("--dec", dest="dec", type=float, default=None,
                       help='DEC of the object (used for querying reference image)')
    parser.add_argument("--filter", dest="filt", default=None,
                       help='Filter of the image')
    parser.add_argument("--sigma_cut", dest="sigma_cut", type=float, default=5,
                       help='Sigma cut for the difference imaging stage (default: 5)')
    parser.add_argument("--show", dest="show", action="store_true",
                       default=False, help='Show results (default: False)')
    parser.add_argument("--download_DES_temp", dest="download_DES_temp", action="store_true",
                       default=False, help='Download DES template (default: False)')
    parser.add_argument("--make_psf_temp", dest="make_psf_temp", action="store_true",
                       default=False, help='Make PSF for template (default: False)')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.download_DES_temp and (args.ra is None or args.dec is None):
        parser.error("--ra and --dec are required when --download_DES_temp is used")
    
    if not args.download_DES_temp and args.tempimg is None:
        parser.error("--tempimg is required when not using --download_DES_temp")
    
    # Get filter from image header if not provided
    if args.filt is None:
        from astropy.io import fits
        filename = os.path.join(args.imgdir, args.targimg)
        try:
            hdul = fits.open(filename)
            header = hdul[0].header
            args.filt = header.get('FILTBAND', 'r')  # Default to 'r' if not found
            hdul.close()
        except Exception as e:
            print(f"Warning: Could not read filter from header: {e}")
            args.filt = 'r'  # Default fallback
    
    # Perform image subtraction
    try:
        output_filename = perform_image_subtraction(
            science_filename=args.targimg,
            template_filename=args.tempimg,
            workdir=args.imgdir,
            ra=args.ra,
            dec=args.dec,
            filt=args.filt,
            sigma_cut=args.sigma_cut,
            show=args.show,
            download_DES_temp=args.download_DES_temp,
            make_psf_temp=args.make_psf_temp
        )
        print(f"Image subtraction completed successfully!")
        print(f"Output file: {output_filename}")
        
    except Exception as e:
        print(f"Error during image subtraction: {e}")
        raise


if __name__ == "__main__":
    main() 