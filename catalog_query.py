from lsst.daf.butler import Butler, Timespan
from lsst.rsp import get_tap_service
import lsst.geom
import lsst.afw.display as afwDisplay
from astropy.time import Time
import numpy as np
import lsst.geom as geom
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord

def gaia3cat(ra, dec, radius_arcmin=10, mag_limit=23, nrows=500):
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


ra = 58.335054  #53.076
dec = -48.750303 #-28.110
band = "r"
time1 = Time("2024-12-01T00:00:00.0", format="isot", scale="tai")
time2 = Time("2024-12-10T00:00:00.0", format="isot", scale="tai")
timespan = Timespan(time1, time2)

dataset_refs = butler.query_datasets("visit_image",
                                     where="band.name = :band AND \
                                     visit.timespan OVERLAPS :timespan AND \
                                     visit_detector_region.region OVERLAPS POINT(:ra, :dec)",
                                     bind={"band": band, "timespan": timespan,
                                           "ra": ra, "dec": dec},
                                     order_by=["visit.timespan.begin"])