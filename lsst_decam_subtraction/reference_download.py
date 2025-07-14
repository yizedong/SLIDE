"""
Reference image downloading utilities for LSST DECam image subtraction.
"""

import numpy as np
import requests
from astropy.table import Table
from astropy.nddata import CCDData
from astropy.wcs import WCS
from astropy.io import fits
from io import BytesIO
import warnings
from pyvo.dal import sia
from numpy.core.defchararray import startswith
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
from astropy import units as u
from .lsst_utils import astropy_world_to_pixel


def download_decam_reference(ra, dec, fov=0.2, filt='g', saveas=None):
    """
    Download DECam reference image from NOIRLab Data Lab.
    
    Parameters
    ----------
    ra, dec : float
        Coordinates in degrees
    fov : float, optional
        Field of view in degrees (default: 0.2)
    filt : str, optional
        Filter band (default: 'g')
    saveas : str, optional
        Path to save file (default: do not save)
        
    Returns
    -------
    CCDData
        Image data with WCS
    """
    DEF_ACCESS_URL = "https://datalab.noirlab.edu/sia/des_dr2"
    svc_des_dr2 = sia.SIAService(DEF_ACCESS_URL)
    imgTable = svc_des_dr2.search((ra,dec), (fov/np.cos(dec*np.pi/180), fov), verbosity=2).to_table()
    sel = (imgTable['proctype'].astype(str)=='Stack') & (imgTable['prodtype'].astype(str)=='image') & (startswith(imgTable['obs_bandpass'].astype(str),filt))
    row = imgTable[sel][0]
    url = row['access_url']
    response = requests.get(url)
    if response.status_code == 200:
        hdul = fits.open(BytesIO(response.content))
        #hdul.info()  # See FITS file structure
        data = hdul[0].data
        header = hdul[0].header
        hdul.close()
        if saveas is not None:
            with fits.open(BytesIO(response.content)) as hdul:
                hdul.info()  # View structure
                hdul.writeto(saveas, overwrite=True)  # Save to disk
    else:
        print(f"Failed to fetch FITS file: {response.status_code}")
    
    ccddata = CCDData(data, wcs=WCS(header), unit='adu')
    ccddata.meta['SATURATE'] = header['SATURATE']
    return ccddata

def gaia3cat_old(ra, dec, ccddata, band='r', radius_arcmin=10, mag_limit=16.5, pm_limit=50, nrows=500):
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

def gaia3cat(ra, dec, ccddata=None, band='r', radius_arcmin=10, mag_limit=21.0, pm_limit=None, nrows=300):
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
    mask = np.isfinite(catalog[mag_column]) & (catalog[mag_column] > mag_limit)
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