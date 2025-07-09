import imp
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


def query_lsst_visits(ra, dec, band, time1, time2):
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
    butler = Butler("dp1", collections="LSSTComCam/DP1")
    timespan = Timespan(time1, time2)
    
    dataset_refs = butler.query_datasets("visit_image",
                                         where="band.name = :band AND \
                                         visit.timespan OVERLAPS :timespan AND \
                                         visit_detector_region.region OVERLAPS POINT(:ra, :dec)",
                                         bind={"band": band, "timespan": timespan,
                                               "ra": ra, "dec": dec},
                                         order_by=["visit.timespan.begin"])
    return dataset_refs

def lsst_visit_to_ccddata(visit_image, saveas=None):
    data = visit_image.image.array
    header = fits.Header(visit_image.getWcs().getFitsMetadata().toDict())
    wcs = WCS(header)
    ccddata = CCDData(data, wcs=wcs, unit='adu')
    if saveas is not None:
        visit_image.writeFits(saveas)
    return ccddata

def lsst_visit_to_psf(visit_image, ra, dec, saveas=None):
    x, y = lsst_world_to_pixel(ra, dec, visit_image)
    psf = visit_image.getPsf()
    xy = lsst.geom.Point2D(x, y)
    psf_image = psf.computeImage(xy)
    if saveas is not None:
        psf_image.writeFits(saveas)
    ccddata = CCDData(psf_image.array, unit='adu')
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

def astropy_world_to_pixel(ra, dec, wcs):
    coord = SkyCoord(ra, dec, unit='deg')
    x, y = skycoord_to_pixel(coord, wcs, origin=0)
    return x, y

def astropy_pixel_to_world(x, y, wcs):
    skycoord = pixel_to_skycoord(x, y, wcs, origin=0)
    return skycoord.ra.deg, skycoord.dec.deg
