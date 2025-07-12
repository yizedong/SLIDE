from lsst.daf.butler import Butler, Timespan
from lsst.rsp import get_tap_service
import lsst.geom
import lsst.afw.display as afwDisplay
from astropy.time import Time
import numpy as np
import lsst.geom as geom
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord
from astropy.nddata import CCDData
from astropy.io import fits
from astropy.wcs import WCS
from photutils.aperture import SkyCircularAnnulus, SkyCircularAperture, ApertureStats
from photutils.background import Background2D
from photutils.psf import ImagePSF
from astropy.table import QTable
from photutils.psf import IterativePSFPhotometry
from photutils.background import LocalBackground, MMMBackground
from photutils.detection import DAOStarFinder
from photutils.psf import PSFPhotometry
from astropy.nddata import Cutout2D


def query_lsst_visits(butler, ra, dec, band, time1, time2):
    """
    Query LSST Butler for visit images within a time range and sky region.
    
    Parameters
    ----------
    butler : lsst.daf.butler.Butler
        Butler instance for data access
    ra : float
        Right ascension in degrees
    dec : float
        Declination in degrees
    band : str
        Filter band name
    time1 : astropy.time.Time
        Start time
    time2 : astropy.time.Time
        End time
        
    Returns
    -------
    list
        Dataset references matching the criteria
    """
    #butler = Butler("dp1", collections="LSSTComCam/DP1")
    time1 = Time(time1, format="isot", scale="tai")
    time2 = Time(time2, format="isot", scale="tai")
    timespan = Timespan(time1, time2)
    
    dataset_refs = butler.query_datasets("visit_image",
                                         where="band.name = :band AND \
                                         visit.timespan OVERLAPS :timespan AND \
                                         visit_detector_region.region OVERLAPS POINT(:ra, :dec)",
                                         bind={"band": band, "timespan": timespan,
                                               "ra": ra, "dec": dec},
                                         order_by=["visit.timespan.begin"])
    return dataset_refs

def lsst_bad_mask(visit_image, mask_type = []):
    mask = visit_image.mask.array
    if len(mask_type) == 0:
        bitmask = visit_image.mask.getPlaneBitMask([
                "BAD",          # detector bad pixel
                "SAT",          # saturated pixel
                "INTRP",        # interpolated pixel
                "CR",           # cosmic ray
                "EDGE",         # edge of the detector
                #"DETECTED",     # source detected
                #"DETECTED_NEGATIVE",  # negative source detected
                "SUSPECT",      # suspicious pixel
                "NO_DATA",      # no data available
                "SENSOR_EDGE",  # sensor edge
                "CLIPPED",      # pixel in clipped area
                "CROSSTALK",    # crosstalk
                #"NOT_DEBLENDED",# not deblended
                "UNMASKEDNAN",  # unmasked NaN
                #"VIGNETTED",    # vignetting
                "STREAK",       # streak from satellite or airplane
            ])
    else:
        bitmask = visit_image.mask.getPlaneBitMask(mask_type)
    image_mask = (visit_image.mask.array & bitmask) != 0
    return image_mask

def lsst_visit_to_ccddata(visit_image, saveas=None):
    data = visit_image.image.array
    header = fits.Header(visit_image.getWcs().getFitsMetadata().toDict())
    wcs = WCS(header)
    mask = lsst_bad_mask(visit_image)
    
    sat_mask = lsst_bad_mask(visit_image, mask_type = ['SAT'])
    masked_values = data[sat_mask] 
    finite_values = masked_values[np.isfinite(masked_values)]
    SATURATE = np.nanmax(finite_values) if finite_values.size > 0 else 60000
    
    ccddata = CCDData(data, wcs=wcs, unit='adu')
    ccddata.mask = mask
    ccddata.meta['SATURATE'] = SATURATE
    if saveas is not None:
        visit_image.writeFits(saveas)
    return ccddata

def get_visit_fwhm(visit_image, ra, dec):
    SIGMA_TO_FWHM = 2.0*np.sqrt(2.0*np.log(2.0))
    x, y = lsst_world_to_pixel(ra, dec, visit_image)
    
    info_calexp = visit_image.getInfo()
    psf_calexp = info_calexp.getPsf()
    
    point_image = lsst.geom.Point2D(x, y)
    psf_shape = psf_calexp.computeShape(point_image)
    psf_sigma = psf_shape.getDeterminantRadius()
    psf_fwhm = psf_sigma * SIGMA_TO_FWHM
    return psf_fwhm

def lsst_visit_to_psf(visit_image, ra, dec):
    x, y = lsst_world_to_pixel(ra, dec, visit_image)
    psf = visit_image.getPsf()
    xy = lsst.geom.Point2D(x, y)
    psf_image = psf.computeImage(xy)
    ccddata = CCDData(psf_image.array, unit='adu')
    return ccddata


def lsst_visit_to_psf_median(visit_image, ra, dec, cutout_size=(2000, 2000), sample_number=20):
    if isinstance(cutout_size, (int, float)):
        cutout_size = (int(cutout_size), int(cutout_size))
    else:
        cutout_size = (int(cutout_size[0]), int(cutout_size[1]))

    data = visit_image.image.array
    ny, nx = data.shape

    # Get WCS and center in pixel coordinates
    header = fits.Header(visit_image.getWcs().getFitsMetadata().toDict())
    wcs = WCS(header)
    xcen, ycen = astropy_world_to_pixel(ra, dec, wcs)

    # Shift center inward if too close to edge
    half_x, half_y = cutout_size[0] // 2, cutout_size[1] // 2
    xcen = np.clip(xcen, half_x, nx - half_x)
    ycen = np.clip(ycen, half_y, ny - half_y)

    # Sample PSFs
    psf = visit_image.getPsf()
    psf_stars = []

    for dx in np.linspace(-cutout_size[0] / 2, cutout_size[0] / 2, sample_number):
        for dy in np.linspace(-cutout_size[1] / 2, cutout_size[1] / 2, sample_number):
            x_sample, y_sample = xcen + dx, ycen + dy
            xy = lsst.geom.Point2D(x_sample, y_sample)
            psf_image = psf.computeImage(xy).getArray()
            psf_stars.append(psf_image)

    if len(psf_stars) == 0:
        raise RuntimeError("No valid PSFs found in the region.")

    # Stack and normalize
    stacked = np.median(np.stack(psf_stars, axis=0), axis=0)
    psf_final = stacked / np.sum(stacked)

    ccddata = CCDData(psf_final, unit='adu')

    return ccddata
        


def lsst_world_to_pixel(ra, dec, visit_image):
    wcs = visit_image.wcs
    if isinstance(ra, float):
        spherePoint = lsst.geom.SpherePoint(ra*geom.degrees, dec*geom.degrees)
        pixel_coord = wcs.skyToPixel(spherePoint)
        x = pixel_coord.getX()
        y = pixel_coord.getY()
    else:
        sphere_points = [geom.SpherePoint(_ra, _dec, geom.degrees) for _ra, _dec in zip(ra, dec)]
        pixel_coords = [wcs.skyToPixel(sp) for sp in sphere_points]
        x = [pt.getX() for pt in pixel_coords]
        y = [pt.getY() for pt in pixel_coords]
    return x, y

def lsst_pixel_to_world(x, y, visit_image):
    wcs = visit_image.wcs
    if isinstance(x, int):
        pixel = geom.Point2D(x, y)
        sky_coord = wcs.pixelToSky(pixel)
        ra = sky_coord.getRa().asDegrees()
        dec = sky_coord.getDec().asDegrees()
    else:
        pixel_coords = [geom.Point2D(px, py) for px, py in zip(x, y)]
        sky_coords = [wcs.pixelToSky(p) for p in pixel_coords]
        ra = [c.getRa().asDegrees() for c in sky_coords]
        dec = [c.getDec().asDegrees() for c in sky_coords]
    return ra, dec

def astropy_world_to_pixel(ra, dec, wcs, origin = 0):
    coord = SkyCoord(ra, dec, unit='deg')
    x, y = skycoord_to_pixel(coord, wcs, origin=origin)
    return x, y

def astropy_pixel_to_world(x, y, wcs, origin=0):
    skycoord = pixel_to_skycoord(x, y, wcs, origin=origin)
    return skycoord.ra.deg, skycoord.dec.deg


def safe_cutout2d(visit_image, ra, dec, cutout_size=(2000, 2000)):
    """
    Create a Cutout2D that contains the original RA/Dec, shifting inward if near the edge.

    Parameters:
        data (2D ndarray): Image array.
        wcs (astropy.wcs.WCS): WCS for the image.
        ra, dec (float): Center coordinates in degrees.
        cutout_size (tuple): Desired cutout size in pixels.

    Returns:
        Cutout2D object.
    """
    if isinstance(cutout_size, (int, float)):
        cutout_size = (int(cutout_size), int(cutout_size))
    else:
        cutout_size = (int(cutout_size[0]), int(cutout_size[1]))

    data = visit_image.image.array
    mask = lsst_bad_mask(visit_image)
    sat_mask = lsst_bad_mask(visit_image, mask_type = ['SAT'])
    masked_values = data[sat_mask] 
    finite_values = masked_values[np.isfinite(masked_values)]
    SATURATE = np.nanmax(finite_values) if finite_values.size > 0 else 60000
    #SATURATE = np.max(data[sat_mask==True])
    header = fits.Header(visit_image.getWcs().getFitsMetadata().toDict())
    wcs = WCS(header)
        
    xcen, ycen = astropy_world_to_pixel(ra, dec, wcs)

    # Half size
    half_x = cutout_size[0] // 2
    half_y = cutout_size[1] // 2
    ny, nx = data.shape

    # Shift center inward if near edges
    xcen = np.clip(xcen, half_x, nx - half_x)
    ycen = np.clip(ycen, half_y, ny - half_y)

    pixel_center = (xcen, ycen)
    cutout = Cutout2D(data, position=pixel_center, size=cutout_size, wcs=wcs, mode='trim')
    cutout_mask = Cutout2D(mask, position=pixel_center, size=cutout_size, wcs=wcs, mode='trim')
    
    ccddata = CCDData(cutout.data, wcs=cutout.wcs, unit='adu')
    ccddata.mask = cutout_mask.data
    ccddata.meta['SATURATE'] = SATURATE
    return ccddata

def lsst_cutout_to_ccddata(cutout, saveas=None):
    ccddata = CCDData(cutout.data, wcs=cutout.wcs, unit='adu')
    return ccddata

def forced_phot(ra, dec, image, wcs, psf_data):
    xx, yy = np.mgrid[:psf_data.shape[0], :psf_data.shape[1]]
    psf_model = ImagePSF(psf_data/np.sum(psf_data), x_0=psf_data.shape[0]/2, y_0=psf_data.shape[1]/2, flux=1)
    bkgstat = MMMBackground()
    localbkg_estimator = LocalBackground(5, 10, bkgstat)
    init_params = QTable()
    _x, _y = astropy_world_to_pixel(ra, dec, wcs)
    init_params['x'] = [_x]
    init_params['y'] = [_y]
    
    fit_shape = (5,5)
    psfphot = PSFPhotometry(psf_model, fit_shape,
                            aperture_radius=5,
                            localbkg_estimator=None)
    phot = psfphot(image, init_params=init_params)
    try:
        flux_njy = phot[0]['flux_fit'].value
        flux_err = phot[0]['flux_err'].value
    except:
        flux_njy = phot[0]['flux_fit']
        flux_err = phot[0]['flux_err']
    mag = -2.5 * np.log10(flux_njy / 3631e9)
    magerr = (2.5 / np.log(10)) * (flux_err / flux_njy)
    upper_limit = -2.5 * np.log10(flux_err * 5 / 3631e9)

    return flux_njy, flux_err, mag, magerr, upper_limit