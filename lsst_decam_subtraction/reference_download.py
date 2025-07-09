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
        hdul.info()  # See FITS file structure
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
    return ccddata

def gaia3cat(ra, dec, radius_arcmin=10, mag_limit=23, nrows=500):
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
        Magnitude limit (default: 23)
    nrows : int, optional
        Maximum number of rows to return (default: 500)
        
    Returns
    -------
    astropy.table.Table
        Gaia catalog results
    """
    coord = SkyCoord(ra=ra, dec=dec, unit='deg')
    query = f"""
    SELECT TOP {nrows} *
    FROM gaiadr3.gaia_source
    WHERE 1=CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', {coord.ra.degree}, {coord.dec.degree}, {radius_arcmin/60.0})
    )
    AND phot_g_mean_mag < {mag_limit}
    """
    job = Gaia.launch_job_async(query)
    return job.get_results()