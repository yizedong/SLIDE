"""
Image processing utilities.
"""

import numpy as np
import matplotlib.pyplot as plt
from reproject import reproject_interp, reproject_adaptive
from astropy.io import fits
from astropy.nddata import CCDData, NDData
from astropy.wcs import WCS, _wcs
from astropy.visualization import PercentileInterval, ImageNormalize
from photutils.psf import EPSFBuilder, extract_stars
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
import scipy.optimize
import warnings
from astropy.wcs.utils import fit_wcs_from_points
import sep
from .lsst_utils import astropy_world_to_pixel

def read_with_datasec(filename, hdu=0):
    """
    Read FITS file with optional DATASEC handling.
    
    Parameters
    ----------
    filename : str
        Path to FITS file
    hdu : int, optional
        HDU index to read (default: 0)
        
    Returns
    -------
    CCDData
        Image data with WCS information
    """
    ccddata = CCDData.read(filename, format='fits', unit='adu', hdu=hdu)
    if 'datasec' in ccddata.meta:
        jmin, jmax, imin, imax = eval(ccddata.meta['datasec'].replace(':', ','))
        ccddata = ccddata[imin-1:imax, jmin-1:jmax]
    return ccddata


def get_ccd_bbox(ccddata):
    """
    Get bounding box coordinates for a CCD image.
    
    Parameters
    ----------
    ccddata : CCDData
        Image data with WCS
        
    Returns
    -------
    tuple
        (ra_min, dec_min, ra_max, dec_max) in degrees
    """
    corners = [[0.], [0.5], [1.]] * np.array(ccddata.shape)[::-1]
    (ra_min, dec_min), (ra_ctr, dec_ctr), (ra_max, dec_max) = ccddata.wcs.all_pix2world(corners, 0.)
    if ra_min > ra_max:
        ra_min, ra_max = ra_max, ra_min
    if dec_min > dec_max:
        dec_min, dec_max = dec_max, dec_min
    max_size_dec = 0.199
    if dec_max - dec_min > max_size_dec:
        dec_min = dec_ctr - max_size_dec / 2.
        dec_max = dec_ctr + max_size_dec / 2.
    return ra_min, dec_min, ra_max, dec_max


def make_psf(data, catalog, show=False, boxsize=25, oversampling=1, center_accuracy=0.001):
    """
    Create PSF model from catalog stars.
    
    Parameters
    ----------
    data : CCDData
        Image data
    catalog : astropy.table.Table
        Star catalog with ra, dec columns
    show : bool, optional
        Show PSF visualization (default: False)
    boxsize : int, optional
        Size of star cutouts (default: 25)
        
    Returns
    -------
    tuple
        (epsf, fitted_stars)
    """
    _catalog = catalog.copy()
    coords = SkyCoord(_catalog['ra'], _catalog['dec'], unit='deg')
    _catalog['x'], _catalog['y'] = skycoord_to_pixel(coords, data.wcs, origin=0)
    
    bkg = np.nanmedian(data.data)
    nddata = NDData(data.data - bkg)

    stars = extract_stars(nddata, _catalog, size=boxsize)
    epsf_builder = EPSFBuilder(oversampling=oversampling, center_accuracy=center_accuracy)
    epsf, fitted_stars = epsf_builder(stars)
    
    if show:
        plt.figure()
        plt.imshow(epsf.data)
        plot_stars(fitted_stars)

    return epsf, fitted_stars

def find_catalog_stars(image, wcs, catalog):
    try:
        data = image.astype(np.float32)
    except:
        data = image.byteswap().newbyteorder()
    bkg = sep.Background(data)
    data_sub = data - bkg

    # 2. Extract all sources
    sources = sep.extract(data_sub, thresh=1.5, err=bkg.globalrms)

    known_x, known_y = astropy_world_to_pixel(catalog['ra'], catalog['dec'], wcs)

    matched_sources = []
    for x0, y0 in zip(known_x, known_y):
        distances = np.hypot(sources['x'] - x0, sources['y'] - y0)
        closest_idx = np.argmin(distances)
        matched_sources.append(sources[closest_idx])

    matched_sources = np.array(matched_sources)
    return matched_sources

def plot_stars(stars):
    """
    Plot extracted stars for visualization.
    
    Parameters
    ----------
    stars : list
        List of star cutouts
    """
    nrows = int(np.ceil(len(stars) ** 0.5))
    fig, axarr = plt.subplots(nrows, nrows, figsize=(20, 20), squeeze=True)
    for ax, star in zip(axarr.ravel(), stars):
        ax.imshow(star)
        ax.plot(star.cutout_center[0], star.cutout_center[1], 'r+')


def update_wcs(wcs, p):
    """
    Update WCS with transformation parameters.
    
    Parameters
    ----------
    wcs : WCS
        WCS object to update
    p : array-like
        Transformation parameters [dx, dy, rotation, scale]
    """
    wcs.wcs.crval += p[:2]
    c, s = np.cos(p[2]), np.sin(p[2])
    if wcs.wcs.has_cd():
        wcs.wcs.cd = wcs.wcs.cd @ np.array([[c, -s], [s, c]]) * p[3]
    if wcs.wcs.has_pc():
        wcs.wcs.pc = wcs.wcs.pc @ np.array([[c, -s], [s, c]]) * p[3]



def wcs_offset(p, radec, xy, origwcs):
    """
    Calculate RMS offset between catalog and detected positions.
    
    Parameters
    ----------
    p : array-like
        Transformation parameters
    radec : array-like
        Catalog coordinates
    xy : array-like
        Detected pixel coordinates
    origwcs : WCS
        Original WCS
        
    Returns
    -------
    float
        RMS offset in pixels
    """
    wcs = origwcs.deepcopy()
    update_wcs(wcs, p)
    sky = SkyCoord(radec[:, 0], radec[:, 1], unit='deg')
    test_xy = np.column_stack(skycoord_to_pixel(sky, wcs))
    rms = (np.sum((test_xy - xy)**2) / len(radec))**0.5
    return rms


def refine_wcs(wcs, stars, catalog, use_sep=False):
    """
    Refine WCS using star positions.
    
    Parameters
    ----------
    wcs : WCS
        WCS to refine
    stars : object
        Detected stars
    catalog : astropy.table.Table
        Reference catalog
    use_sep : bool, optional
        Use SEP detection format (default: False)
    """
    if use_sep:
        xy = np.array([[star['x'], star['y']] for star in stars])
        t_match = catalog[stars['i']]
    else:
        xy = np.array([star.center for star in stars.all_good_stars])
        t_match = catalog[[star.id_label - 1 for star in stars.all_good_stars]]
    radec = np.array([t_match['raMean'], t_match['decMean']]).T

    res = scipy.optimize.minimize(wcs_offset, [0., 0., 0., 1.], args=(radec, xy, wcs),
                                  bounds=[(-0.01, 0.01), (-0.01, 0.01), (-0.1, 0.1), (0.9, 1.1)])

    orig_rms = wcs_offset([0., 0., 0., 1.], radec, xy, wcs)
    print(' orig_fun: {}'.format(orig_rms))
    print(res)
    update_wcs(wcs, res.x)

def refine_wcs_astropy(image, wcs, catalog, projection='TAN'):
    sky_coords = SkyCoord(ra=catalog['ra'], dec=catalog['dec'])
    matched_sources = find_catalog_stars(image, wcs, catalog)

    # matched pixel positions and sky positions
    new_wcs = fit_wcs_from_points(xy=(matched_sources['x'], matched_sources['y']), world_coords=sky_coords, projection=projection,
                                 sip_degree=5)
    return new_wcs
    

def assemble_reference(refdatas, wcs, shape, ref_global_bkg=0, order='bicubic'):
    """
    Reproject and stack reference images to match science image.
    
    Parameters
    ----------
    refdatas : list
        List of reference CCDData objects
    wcs : WCS
        Target WCS
    shape : tuple
        Target image shape
    ref_global_bkg : float, optional
        Background value for masked pixels (default: 0)
    order : str, optional
        Interpolation order (default: 'bicubic')
        
    Returns
    -------
    CCDData
        Assembled reference image
    """
    refdatas_reprojected = []
    refdata_foot = np.zeros(shape, float)
    for data in refdatas:
        #reprojected, foot = reproject_interp((data.data, data.wcs), wcs, shape, order=order)
        reprojected, foot = reproject_adaptive((data.data, data.wcs), wcs, shape, conserve_flux=True)
        
        refdatas_reprojected.append(reprojected)
        refdata_foot += foot

    refdata_reproj = np.nanmean(refdatas_reprojected, axis=0)
    if np.all(np.isnan(refdata_reproj)):
        raise ValueError("All reprojected reference data is NaN.")
    refdata_reproj[np.isnan(refdata_reproj)] = ref_global_bkg
    refdata = CCDData(refdata_reproj, wcs=wcs, mask=refdata_foot == 0., unit='adu')
    return refdata 