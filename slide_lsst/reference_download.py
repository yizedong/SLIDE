"""
Reference image downloading utilities for LSST DECam image subtraction.
"""

import gzip
import numpy as np
import requests
from astropy.table import Table
from astropy.nddata import CCDData, Cutout2D
from astropy.wcs import WCS
from astropy.io import fits
from io import BytesIO
import warnings
from pyvo.dal import sia, TAPService
from numpy.core.defchararray import startswith
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
from astropy import units as u
from .lsst_utils import astropy_world_to_pixel
from astropy.utils.exceptions import AstropyUserWarning

_NERSC_BASE = "https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/coadd"
_BRICKS_URL = "https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/survey-bricks-dr10-south.fits.gz"
_bricks_cache = None


def _load_dr10_bricks():
    global _bricks_cache
    if _bricks_cache is None:
        response = requests.get(_BRICKS_URL, timeout=120)
        response.raise_for_status()
        with fits.open(BytesIO(gzip.decompress(response.content))) as hdul:
            _bricks_cache = hdul[1].data
    return _bricks_cache


def _find_dr10_brickname(ra, dec):
    bricks = _load_dr10_bricks()
    mask = (
        (bricks['RA1'] <= ra)  & (ra  < bricks['RA2']) &
        (bricks['DEC1'] <= dec) & (dec < bricks['DEC2'])
    )
    hits = bricks[mask]
    if len(hits) == 0:
        raise ValueError(
            f'RA={ra}, Dec={dec} not found in any DR10-South brick. '
            'Check that the position is within the DECam footprint (Dec < 32.375°).'
        )
    return hits['BRICKNAME'][0].strip()


def download_des_reference(ra, dec, fov=0.2, filt='g', saveas=None):
    """
    Download DECam reference image and mask from NOIRLab Data Lab and optionally save to a single FITS file.

    Parameters
    ----------
    ra, dec : float
        Coordinates in degrees.
    fov : float, optional
        Field of view in degrees (default: 0.2).
    filt : str, optional
        Filter band (default: 'g').
    saveas : str, optional
        Path to save combined FITS file (image + mask in HDUs).

    Returns
    -------
    CCDData
        Image data with mask and WCS.
    """
    DEF_ACCESS_URL = "https://datalab.noirlab.edu/sia/des_dr2"
    svc_des_dr2 = sia.SIAService(DEF_ACCESS_URL)
    imgTable = svc_des_dr2.search((ra, dec), (fov / np.cos(np.radians(dec)), fov), verbosity=2).to_table()
    if len(imgTable) == 0:
        warnings.warn("No image entries found for given coordinates in DES.", AstropyUserWarning)
        return None
        
    # Filter for proper image entries
    sel = (
        (imgTable['proctype'].astype(str) == 'Stack') &
        (imgTable['prodtype'].astype(str) == 'image') &
        (startswith(imgTable['obs_bandpass'].astype(str), filt))
    )
    selected = imgTable[sel]
    if len(selected) < 2:
        raise ValueError("Expected at least 2 rows (image + mask), but found fewer.")

    # --- Science image ---
    img_url = selected[0]['access_url']
    response = requests.get(img_url)
    response.raise_for_status()
    with fits.open(BytesIO(response.content)) as hdul:
        data = hdul[0].data
        header = hdul[0].header

    # --- Mask ---
    mask_url = selected[1]['access_url']
    response = requests.get(mask_url)
    response.raise_for_status()
    with fits.open(BytesIO(response.content)) as hdul:
        mask_data = hdul[0].data.astype(bool)

    # --- Save combined FITS file if requested ---
    if saveas:
        hdu_image = fits.PrimaryHDU(data=data, header=header)
        hdu_mask = fits.ImageHDU(data=mask_data.astype(np.uint8), name='mask')  # Save mask as 0/1 integers
        fits.HDUList([hdu_image, hdu_mask]).writeto(saveas, overwrite=True)

    # --- Create CCDData ---
    ccddata = CCDData(data, wcs=WCS(header), unit='adu', mask=mask_data)
    ccddata.meta['SATURATE'] = header.get('SATURATE', 65535)
    ccddata.meta['ref_zp'] = 30 #mag ZP
    return ccddata


def download_decals_reference(ra, dec, fov=0.2, filt='g', saveas=None):
    """
    Download DECam reference image from ls dr9.

    Parameters
    ----------
    ra, dec : float
        Coordinates in degrees.
    fov : float, optional
        Field of view in degrees (default: 0.2).
    filt : str, optional
        Filter band (default: 'g').
    saveas : str, optional
        Path to save combined FITS file (image + mask in HDUs).

    Returns
    -------
    CCDData
        Image data with mask and WCS.
    """
    DEF_ACCESS_URL = "https://datalab.noirlab.edu/sia/coadd/ls_dr9"
    svc_des_dr2 = sia.SIAService(DEF_ACCESS_URL)
    imgTable = svc_des_dr2.search((ra, dec), (fov / np.cos(np.radians(dec)), fov), verbosity=2).to_table()

    # Filter for proper image entries
    sel = (
        (imgTable['proctype'].astype(str) == 'Stack') &
        (imgTable['prodtype'].astype(str) == 'image') &
        (startswith(imgTable['obs_bandpass'].astype(str), filt))
    )
    selected = imgTable[sel]
    #if len(selected) < 2:
    #    raise ValueError("Expected at least 2 rows (image + mask), but found fewer.")

    # --- Science image ---
    img_url = selected[0]['access_url']
    response = requests.get(img_url)
    response.raise_for_status()
    with fits.open(BytesIO(response.content)) as hdul:
        data = hdul[0].data
        header = hdul[0].header

    # --- Mask ---
    mask_data = (data == 0).astype(bool)

    # --- Save combined FITS file if requested ---
    if saveas:
        hdu_image = fits.PrimaryHDU(data=data, header=header)
        hdu_mask = fits.ImageHDU(data=mask_data.astype(np.uint8), name='mask')  # Save mask as 0/1 integers
        fits.HDUList([hdu_image, hdu_mask]).writeto(saveas, overwrite=True)

    # --- Create CCDData ---
    ccddata = CCDData(data, wcs=WCS(header), unit='adu', mask=mask_data)
    ccddata.meta['SATURATE'] = header.get('SATURATE', 65535) #not sure how to get the correct saturation level for ld dr9
    ccddata.meta['ref_zp'] = header.get('MAGZERO') #mag ZP
    return ccddata

def download_decals_dr10_reference(ra, dec, fov=0.2, filt='g', saveas=None):
    """
    Download DECam reference image from ls dr9.

    Parameters
    ----------
    ra, dec : float
        Coordinates in degrees.
    fov : float, optional
        Field of view in degrees (default: 0.2).
    filt : str, optional
        Filter band (default: 'g').
    saveas : str, optional
        Path to save combined FITS file (image + mask in HDUs).

    Returns
    -------
    CCDData
        Image data with mask and WCS.
    """
    DEF_ACCESS_URL = "https://datalab.noirlab.edu/sia/ls_dr10"
    svc_des_dr2 = sia.SIAService(DEF_ACCESS_URL)
    imgTable = svc_des_dr2.search((ra, dec), (fov / np.cos(np.radians(dec)), fov), verbosity=2).to_table()

    # Filter for proper image entries
    sel = (
        (imgTable['proctype'].astype(str) == 'Stack') &
        (imgTable['prodtype'].astype(str) == 'image') &
        (startswith(imgTable['obs_bandpass'].astype(str), filt))
    )
    selected = imgTable[sel]
    #if len(selected) < 2:
    #    raise ValueError("Expected at least 2 rows (image + mask), but found fewer.")

    # --- Science image ---
    img_url = selected[0]['access_url']
    response = requests.get(img_url)
    response.raise_for_status()
    with fits.open(BytesIO(response.content)) as hdul:
        data = hdul[0].data
        header = hdul[0].header

    # --- Mask ---
    mask_data = (data == 0).astype(bool)

    # --- Save combined FITS file if requested ---
    if saveas:
        hdu_image = fits.PrimaryHDU(data=data, header=header)
        hdu_mask = fits.ImageHDU(data=mask_data.astype(np.uint8), name='mask')  # Save mask as 0/1 integers
        fits.HDUList([hdu_image, hdu_mask]).writeto(saveas, overwrite=True)

    # --- Create CCDData ---
    ccddata = CCDData(data, wcs=WCS(header), unit='adu', mask=mask_data)
    ccddata.meta['SATURATE'] = header.get('SATURATE', 65535) #not sure how to get the correct saturation level for ld dr9
    ccddata.meta['ref_zp'] = header.get('MAGZERO') #mag ZP
    return ccddata
    

def gaia3cat_old_1(ra, dec, ccddata, band='r', radius_arcmin=10, mag_limit=16.5, pm_limit=50, nrows=500):
    """
    Query Gaia DR3 catalog for stars within a specified radius.
    
    Parameters
    ----------
    ra : float
        Right ascension in degrees
    dec : float
        Declination in degrees
    radius_arcmin : float, optional
        Search radius in arcminutes (default: 10)
    mag_limit : float, optional
        Magnitude limit (default: 16.5)
    nrows : int, optional
        Maximum number of rows to return (default: 500)
        
    Returns
    -------
    astropy.table.Table
        Gaia catalog results
    """
    if band in ['u', 'g']:
        gaia_filt = 'bp'
    elif band in ['r', 'i', 'z']:
        gaia_filt = 'rp'
    else:
        gaia_filt = 'g' 
    coord = SkyCoord(ra=ra, dec=dec, unit='deg')
    query = f"""
    SELECT TOP {nrows} *
    FROM gaiadr3.gaia_source
    WHERE 1=CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', {coord.ra.degree}, {coord.dec.degree}, {radius_arcmin/60.0})
    )
    AND phot_{gaia_filt}_mean_mag > {mag_limit}
    AND SQRT(POWER(pmra, 2) + POWER(pmdec, 2)) < {pm_limit}
    """
    job = Gaia.launch_job_async(query)
    catalog = job.get_results()
    
    # Filter stars within image bounds
    nx, ny = ccddata.shape
    x, y = astropy_world_to_pixel(catalog['ra'], catalog['dec'], ccddata.wcs)
    margin = 25 // 2
    
    mask = (
        (x > margin) & (x < nx - margin) &
        (y > margin) & (y < ny - margin)
    )
    catalog = catalog[mask]
    return catalog

def gaia3cat_old(ra, dec, ccddata=None, band='r', radius_arcmin=10, mag1=17, mag2=22, pm_limit=None, nrows=300):
    """
    Query Gaia DR3 catalog using a rectangular (box) query for stars around a given sky position.

    Parameters
    ----------
    ra : float
        Right ascension in degrees.
    dec : float
        Declination in degrees.
    ccddata : CCDData, optional
        If provided, stars outside the CCD frame will be excluded.
    band : str, optional
        Observing band to choose Gaia magnitude column ('bp', 'rp', or 'g'). Default: 'r' -> 'rp'.
    radius_arcmin : float, optional
        Search box size in arcminutes (default: 10). Width is declination-corrected.
    mag_limit : float, optional
        Maximum magnitude to include (default: 21.0).
    pm_limit : float or None, optional
        If given, exclude sources with total proper motion above this value (mas/yr).
    nrows : int, optional
        Maximum number of rows to retrieve from Gaia (default: 300).

    Returns
    -------
    catalog : astropy.table.Table
        Filtered Gaia catalog.
    """
    Gaia.ROW_LIMIT = nrows

    coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    height = u.Quantity(radius_arcmin*2, u.arcmin)
    width = u.Quantity(radius_arcmin*2 / np.cos(np.radians(dec)), u.arcmin)

    catalog = Gaia.query_object_async(coordinate=coord, width=width, height=height)

    # Choose Gaia magnitude column based on band
    if band in ['u', 'g']:
        mag_column = 'phot_bp_mean_mag'
    elif band in ['r', 'i', 'z']:
        mag_column = 'phot_rp_mean_mag'
    else:
        mag_column = 'phot_g_mean_mag'

    # Apply filters
    mask = np.isfinite(catalog[mag_column]) & (catalog[mag_column] >= mag1) & (catalog[mag_column] <= mag2)
    mask &= catalog['astrometric_excess_noise_sig'] < 2
    if pm_limit is not None:
        total_pm = np.hypot(catalog['pmra'], catalog['pmdec'])
        mask &= total_pm < pm_limit

    catalog = catalog[mask]

    # Clip to CCD footprint if ccddata is provided
    if ccddata is not None:
        nx, ny = ccddata.shape
        x, y = astropy_world_to_pixel(catalog['ra'], catalog['dec'], ccddata.wcs)
        margin = 25 // 2
        spatial_mask = (x > margin) & (x < nx - margin) & (y > margin) & (y < ny - margin)
        catalog = catalog[spatial_mask]

    return catalog

def gaia3cat(ra, dec, ccddata=None, band='r', radius_arcmin=10, mag1=17, mag2=22, pm_limit=None, nrows=300):
    """
    Query Gaia DR3 from the NOIRLab Data Lab TAP service instead of the ESA archive.

    Parameters
    ----------
    ra : float
        Right ascension in degrees.
    dec : float
        Declination in degrees.
    ccddata : CCDData, optional
        If provided, stars outside the CCD frame will be excluded.
    band : str, optional
        Observing band to choose Gaia magnitude column ('bp', 'rp', or 'g'). Default: 'r' -> 'rp'.
    radius_arcmin : float, optional
        Search box half-size in arcminutes (default: 10). Width is declination-corrected.
    mag1, mag2 : float, optional
        Magnitude range to include (default: 17 to 22).
    pm_limit : float or None, optional
        If given, exclude sources with total proper motion above this value (mas/yr).
    nrows : int, optional
        Maximum number of rows to keep, closest to (ra, dec) first (default: 300).

    Returns
    -------
    catalog : astropy.table.Table
        Filtered Gaia catalog, with 'ra' and 'dec' columns.
    """
    # Choose Gaia magnitude column based on band
    if band in ['u', 'g']:
        mag_column = 'phot_bp_mean_mag'
    elif band in ['r', 'i', 'z']:
        mag_column = 'phot_rp_mean_mag'
    else:
        mag_column = 'phot_g_mean_mag'

    half_dec = radius_arcmin / 60.0
    half_ra = half_dec / np.cos(np.radians(dec))

    # Data Lab's TAP does not take the ADQL geometry functions, so cut on ra/dec
    # directly and split the range when the box straddles ra = 0.
    dec1, dec2 = dec - half_dec, dec + half_dec
    ra1, ra2 = ra - half_ra, ra + half_ra
    if ra1 < 0 or ra2 > 360:
        ra_cut = f'(ra >= {ra1 % 360} OR ra <= {ra2 % 360})'
    else:
        ra_cut = f'ra BETWEEN {ra1} AND {ra2}'

    columns = ('ra, dec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, '
               'pmra, pmdec, astrometric_excess_noise_sig')
    query = f"""
    SELECT {columns}
    FROM gaia_dr3.gaia_source
    WHERE {ra_cut}
      AND dec BETWEEN {dec1} AND {dec2}
      AND {mag_column} BETWEEN {mag1} AND {mag2}
      AND astrometric_excess_noise_sig < 2
    """
    catalog = TAPService('https://datalab.noirlab.edu/tap').search(query).to_table()

    if len(catalog) == 0:
        warnings.warn('No Gaia DR3 sources found for the given box.', AstropyUserWarning)
        return catalog

    if pm_limit is not None:
        total_pm = np.hypot(catalog['pmra'], catalog['pmdec'])
        catalog = catalog[total_pm < pm_limit]

    # Data Lab cannot sort by distance in the query, so keep the nearest nrows here
    if len(catalog) > nrows:
        separation = SkyCoord(ra, dec, unit='deg').separation(
            SkyCoord(catalog['ra'], catalog['dec'], unit='deg'))
        catalog = catalog[np.argsort(separation)[:nrows]]

    # Clip to CCD footprint if ccddata is provided
    if ccddata is not None:
        nx, ny = ccddata.shape
        x, y = astropy_world_to_pixel(catalog['ra'], catalog['dec'], ccddata.wcs)
        margin = 25 // 2
        spatial_mask = (x > margin) & (x < nx - margin) & (y > margin) & (y < ny - margin)
        catalog = catalog[spatial_mask]

    return catalog

def des2cat(ra, dec, ccddata, band='r', radius_arcmin=11, mag1=16, mag2=21):
    radius = radius_arcmin
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree))
    Vizier.ROW_LIMIT=-1
    Vizier.columns = ['RA_ICRS', 'DE_ICRS', 'Tile',
                      'CoadID', 
                      f'S/G{band}',
                      f'{band}mag0',
                       f'{band}FWHM']
    Vizier.column_filters={'S/Gr': '>0.5',
                            f'{band}mag0': f'>{mag1:.2f} && <={mag2:.2f}'}
    t = Vizier.query_region(coord,
                                 radius=float('{:f}'.format(radius))*u.arcmin,
                                 catalog='II/371')
    des_table = t[0]
    des_table['ra'] = des_table['RA_ICRS']
    des_table['dec'] = des_table['DE_ICRS']

    catalog = des_table
    # Filter stars within image bounds
    nx, ny = ccddata.shape
    x, y = astropy_world_to_pixel(catalog['ra'], catalog['dec'], ccddata.wcs)
    margin = 25 // 2
    
    mask = (
        (x > margin) & (x < nx - margin) &
        (y > margin) & (y < ny - margin)
    )
    catalog = catalog[mask]
    return catalog

def apass2cat(ra, dec, ccddata, band='r', radius_arcmin=11):
    """
    Query APASS DR9 catalog and return sources within image bounds and mag range.
    
    Parameters:
        ra, dec (float): Target sky coordinates (degrees).
        ccddata (CCDData): Image with WCS info.
        band (str): Photometric band, one of 'B', 'V', 'g', 'r', 'i'.
        radius_arcmin (float): Search radius in arcminutes.
        mag1, mag2 (float): Magnitude range.

    Returns:
        Table: Filtered APASS catalog entries within image and mag range.
    """
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
    Vizier.ROW_LIMIT = -1
    
    # Set up Vizier filters
    mag_col = f"{band}mag"
    Vizier.columns = ['RAJ2000', 'DEJ2000', "g'mag", "r'mag", "i'mag"] #, "g'mag", "r'mag", "i'mag"
    #Vizier.column_filters = {mag_col: f">{mag1:.2f} && <={mag2:.2f}"}

    # Query APASS DR9
    result = Vizier.query_region(coord, radius=radius_arcmin*u.arcmin, catalog="II/336/apass9")
    if len(result) == 0:
        return None

    apass_table = result[0]
    apass_table['ra'] = apass_table['RAJ2000']
    apass_table['dec'] = apass_table['DEJ2000']

    # Project to pixel coordinates and filter by image bounds
    nx, ny = ccddata.shape
    x, y = astropy_world_to_pixel(apass_table['ra'], apass_table['dec'], ccddata.wcs)
    margin = 25 // 2
    mask = (x > margin) & (x < nx - margin) & (y > margin) & (y < ny - margin)
    return apass_table[mask]