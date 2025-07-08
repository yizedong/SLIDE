from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from reproject import reproject_interp
import os
from astropy.wcs import WCS, _wcs
from astropy.visualization import PercentileInterval, ImageNormalize
import requests
from astropy.table import Table
from astropy.nddata import CCDData, NDData
#from photutils import psf, EPSFBuilder
from photutils.psf import EPSFBuilder, extract_stars
from io import BytesIO
from PyZOGY.subtract import calculate_difference_image, calculate_difference_image_zero_point, \
    normalize_difference_image, save_difference_image_to_file, calculate_difference_psf, save_difference_psf_to_file
from PyZOGY.image_class import ImageClass
import scipy
import warnings
import sep
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.vizier import Vizier
import warnings
from argparse import ArgumentParser
from pyvo.dal import sia
from astropy.utils.data import download_file
from astropy.utils.data import clear_download_cache
from numpy.core.defchararray import startswith
from astropy.wcs.utils import skycoord_to_pixel
from astroquery.gaia import Gaia


def read_with_datasec(filename, hdu=0):
    ccddata = CCDData.read(filename, format='fits', unit='adu', hdu=hdu)
    if 'datasec' in ccddata.meta:
        jmin, jmax, imin, imax = eval(ccddata.meta['datasec'].replace(':', ','))
        ccddata = ccddata[imin-1:imax, jmin-1:jmax]
    return ccddata


def get_ccd_bbox(ccddata):
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


def get_ps1_catalog(ra_min, dec_min, ra_max, dec_max, mag_max=21., mag_min=16., mag_filter='r'):
    res = requests.get('http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx',
                       params={'cat': 'PS1V3OBJECTS', 'format': 'csv', 'mindet': 25,
                               'bbox': '{},{},{},{}'.format(ra_min, dec_min, ra_max, dec_max)})
    t = Table.read(res.text, format='csv', header_start=1, data_start=2)
    psfmag_key = mag_filter + 'MeanPSFMag'
    is_point_source = t[psfmag_key] - t[mag_filter + 'MeanKronMag'] < 0.05
    mag_cut = (t[psfmag_key] < mag_max) & (t[psfmag_key] > mag_min)
    t_stars = t[is_point_source & mag_cut]
    return t_stars


def make_psf(data, catalog, show=False, boxsize=25):
    catalog = catalog.copy()
    #catalog['x'], catalog['y'] = data.wcs.all_world2pix(catalog['raMean'], catalog['decMean'], 0)
    coords = SkyCoord(catalog['raMean'], catalog['decMean'], unit='deg')
    catalog['x'], catalog['y'] = skycoord_to_pixel(coords, data.wcs, origin=0)
    
    bkg = np.nanmedian(data.data)
    nddata = NDData(data.data - bkg)

    stars = extract_stars(nddata, catalog, size=boxsize)
    epsf_builder = EPSFBuilder(oversampling=1)
    epsf, fitted_stars = epsf_builder(stars)
    
    if show:
        plt.figure()
        plt.imshow(epsf.data)
        plot_stars(fitted_stars)

    return epsf, fitted_stars


def plot_stars(stars):
    nrows = int(np.ceil(len(stars) ** 0.5))
    fig, axarr = plt.subplots(nrows, nrows, figsize=(20, 20), squeeze=True)
    for ax, star in zip(axarr.ravel(), stars):
        ax.imshow(star)
        ax.plot(star.cutout_center[0], star.cutout_center[1], 'r+')


def update_wcs(wcs, p):
    wcs.wcs.crval += p[:2]
    c, s = np.cos(p[2]), np.sin(p[2])
    if wcs.wcs.has_cd():
        wcs.wcs.cd = wcs.wcs.cd @ np.array([[c, -s], [s, c]]) * p[3]
    if wcs.wcs.has_pc():
        wcs.wcs.pc = wcs.wcs.pc @ np.array([[c, -s], [s, c]]) * p[3]


def wcs_offset(p, radec, xy, origwcs):
    wcs = origwcs.deepcopy()
    update_wcs(wcs, p)
    #test_xy = wcs.all_world2pix(radec, 0)
    sky = SkyCoord(radec[:, 0], radec[:, 1], unit='deg')
    test_xy = np.column_stack(skycoord_to_pixel(sky, wcs))
    rms = (np.sum((test_xy - xy)**2) / len(radec))**0.5
    return rms


def refine_wcs(wcs, stars, catalog, use_sep=False):
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


def get_ps1_filename(ra, dec, filt):
    """
    Download Image from PS1 and correct luptitudes back to a linear scale.

    Parameters
    ---------------
    ra, dec : Coordinates in degrees
    filt    : Filter color 'g', 'r', 'i', 'z', or 'y'

    Output
    ---------------
    filename : PS1 image filename
    """

    # Query a center RA and DEC from PS1 in a specified color
    res = requests.get('http://ps1images.stsci.edu/cgi-bin/ps1filenames.py',
                 params={'ra': ra, 'dec': dec, 'filters': filt})
    t = Table.read(res.text, format='ascii')

    return t['filename'][0]


def download_ps1_image(filename, saveas=None):
    """
    Download image from PS1 and correct luptitudes back to a linear scale.

    Parameters
    ---------------
    filename : PS1 image filename (from `get_ps1_filename`)
    saveas   : Path to save template file (default: do not save)

    Output
    ---------------
    ccddata : CCDData format of data with WCS
    """
    res = requests.get('http://ps1images.stsci.edu' + filename)
    hdulist = fits.open(BytesIO(res.content))

    # Linearize from luptitudes
    boffset = hdulist[1].header['boffset']
    bsoften = hdulist[1].header['bsoften']
    data_linear = boffset + bsoften * 2 * np.sinh(hdulist[1].data * np.log(10.) / 2.5)
    warnings.simplefilter('ignore')  # ignore warnings from nonstandard PS1 header keywords
    ccddata = CCDData(data_linear, wcs=WCS(hdulist[1].header), unit='adu')

    # Save the template to file
    if saveas is not None:
        ccddata.write(saveas, overwrite=True)

    return ccddata

def download_decam_reference(ra, dec, fov=0.2, filt='g', saveas=None):
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

    # Save the template to file
    if saveas is not None:
        ccddata.write(saveas, overwrite=True)

    return ccddata

def download_references(ra_min, dec_min, ra_max, dec_max, mag_filter, template_basename=None, catalog=None):
    """
    Download 1 to 4 references from PS1 as necessary to cover full RA & dec range

    Parameters
    ---------------
    ra_min, ra_max   : Minimum and Maximum RA and DEC
    dec_min, dec_max   in units of degrees
    mag_filter       : Filter color 'g', 'r', 'i', 'z', or 'y'
    template_basename: Filename of the output(s), to be suffixed by 0.fits, 1.fits, ...
    catalog          : Catalog to which to align the reference image WCS (default: do not align)

    Output
    ---------------
    refdatas   : List of CCDData objects containing the reference images

    """

    filename0 = get_ps1_filename(ra_min, dec_min, mag_filter)
    filename1 = get_ps1_filename(ra_max, dec_max, mag_filter)
    filename2 = get_ps1_filename(ra_min, dec_max, mag_filter)
    filename3 = get_ps1_filename(ra_max, dec_min, mag_filter)

    filenames = {filename0, filename1, filename2, filename3}
    refdatas = []
    for i, fn in enumerate(filenames):
        if template_basename is not None:
            saveas = template_basename + '{:d}.fits'.format(i)
            print('downloading', saveas)
        else:
            saveas = None
            print('downloading', fn)
        refdata = download_ps1_image(fn, saveas)
        if catalog is not None:
            _, stars = make_psf(refdata, catalog)
            try:
                refine_wcs(refdata.wcs, stars, catalog)
            except _wcs.InvalidTransformError:
                print('WARNING: unable to refine wcs')
            except ValueError:
                print('WARNING: no stars found, unable to refine wcs')
        refdatas.append(refdata)

    return refdatas




def assemble_reference(refdatas, wcs, shape, ref_global_bkg=0, order='bicubic'):
    """Reproject and stack the reference images to match the science image"""
    refdatas_reprojected = []
    refdata_foot = np.zeros(shape, float)
    for data in refdatas:
        reprojected, foot = reproject_interp((data.data, data.wcs), wcs, shape, order=order)
        refdatas_reprojected.append(reprojected)
        refdata_foot += foot

    refdata_reproj = np.nanmean(refdatas_reprojected, axis=0)
    #refdata_reproj[np.isnan(refdata_reproj)] = 0.
    if np.all(np.isnan(refdata_reproj)):
        raise ValueError("All reprojected reference data is NaN.")
    refdata_reproj[np.isnan(refdata_reproj)] = ref_global_bkg #np.nanmedian(refdata_reproj)
    refdata = CCDData(refdata_reproj, wcs=wcs, mask=refdata_foot == 0., unit='adu')
    return refdata

def gaia2file(ra, dec, size=26., mag_limit=18., output='gaia.cat'):

    from astroquery.gaia import Gaia

    warnings.simplefilter('ignore')  # suppress a lot of astroquery warnings

    coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree))
    height = u.Quantity(size, u.arcminute)
    width  = u.Quantity(size/np.cos(dec*np.pi/180.), u.arcminute)
    response = Gaia.query_object_async(coordinate=coord, width=width, height=height)

    response = response[
            (response['phot_g_mean_mag'] < mag_limit) &
            (response['astrometric_excess_noise_sig'] < 2) # filter objects with bad astrometry (e.g. HII regions in galaxies)
    ]
    response['ra'].format ='%16.12f'
    response['dec'].format = '%16.12f'
    response['phot_g_mean_mag'].format = '%.2f'

    gaia_cat = response['ra', 'dec', 'SOURCE_ID', 'phot_g_mean_mag']
    #gaia_cat.write(output, format='ascii.commented_header',
   #         delimiter=' ', overwrite=True)
    return gaia_cat

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

if __name__ == '__main__':
    description = ">> make different images for LSST"
    usage = "%(prog)s [options] "
    parser = ArgumentParser(usage=usage, description=description)
    parser.add_argument("--targimg", dest="targimg", default='None',
            help='the name of the science image')
    parser.add_argument("--tempimg", dest="tempimg", default='None',
            help='the name of the template image')
    parser.add_argument("--imgdir", dest="imgdir", default='./',
            help='the path of the images')
    parser.add_argument("--ra", dest="ra", default=0, type = float,
            help='ra of the object (used for querying the reference image)')
    parser.add_argument("--dec", dest="dec", default=0, type = float,
            help='dec of the object (used for querying the reference image)')
    parser.add_argument("--filter", dest="filt", default=0,
            help='filter of the image')
    parser.add_argument("--sigma_cut", dest="sigma_cut", default=5,
            help='sigma cut for the difference imaing stage')
    parser.add_argument("--show", dest="show", action="store_true",
                      default=False, help=' show result  \t\t\t [%default]')
    parser.add_argument("--download_DES_temp", dest="download_DES_temp", action="store_true",
                      default=False, help=' make the psf using photutils  \t\t\t [%default]')
    parser.add_argument("--make_psf_temp", dest="make_psf_temp", action="store_true",
                      default=False, help=' make the psf using photutils  \t\t\t [%default]')


    args = parser.parse_args()

    #show = False
    workdir = args.imgdir #'./img_sub/'
    template_filename = args.tempimg
    science_filename = args.targimg
    show = args.show
    make_psf_temp = args.make_psf_temp
    sigma_cut = args.sigma_cut
    SN_ra = args.ra
    SN_dec = args.dec
    download_DES_temp = args.download_DES_temp
    filt = args.filt

    # # Read the science image
    filename = os.path.join(workdir, science_filename)
    scidata0 = read_with_datasec(filename)
    nx, ny = scidata0.shape
    half_ra, half_dec = scidata0.wcs.all_pix2world(nx//2, ny//2, 0)
    #catalog = Table.read('gaia.cat', format='ascii')
    print(np.round(half_ra, 3), np.round(half_dec, 3))
    #catalog = gaia2file(ra=np.round(half_ra, 3), dec=np.round(half_dec, 3), size=20., mag_limit=23., output='gaia.cat')
    catalog = gaia3cat(ra=np.round(half_ra, 3), dec=np.round(half_dec, 3), radius_arcmin=10)
    print('gaia catalog size:', len(catalog))
    catalog['raMean'], catalog['decMean'] = catalog['ra'], catalog['dec']
    coords = SkyCoord(catalog['raMean'], catalog['decMean'], unit='deg')
    x, y = skycoord_to_pixel(coords, scidata0.wcs, origin=0)
    #x, y = scidata0.wcs.all_world2pix(catalog['raMean'], catalog['decMean'], 0)
    margin = 25 // 2
    
    mask = (
        (x > margin) & (x < nx - margin) &
        (y > margin) & (y < ny - margin)
    )
    catalog = catalog[mask]

    #fig, ax2 = plt.subplots(1, 1, figsize=(15, 15))
    #norm = ImageNormalize(scidata0.data, PercentileInterval(99.))
    #ax2.imshow(scidata0.data, norm=norm)
    #ax2.plot(x, y, marker='o', mec='r', mfc='none', ls='none')
    #plt.show()
        

    _, sci_stars = make_psf(scidata0, catalog, show=show, boxsize=25)
    scidata = scidata0.copy()
    refine_wcs(scidata.wcs, sci_stars, catalog)
    print(len(catalog))

    filename = os.path.join(workdir, science_filename.replace('.fits', '.psf.fits'))
    sci_psf = read_with_datasec(filename)

    # # Download the reference image
    if download_DES_temp:
        _refdatas = download_decam_reference(ra=SN_ra, dec=SN_dec, fov=0.5, filt=filt)
        print(f'DS template downloaded for RA:{SN_ra} DEC:{SN_dec}')
    else:
        filename = os.path.join(workdir, template_filename)
        _refdatas = read_with_datasec(filename)
    nx, ny = _refdatas.shape
    #half_ra, half_dec = _refdatas.wcs.all_pix2world(nx//2, ny//2, 0)
    half_ra, half_dec = SN_ra, SN_dec
    #catalog = Table.read('gaia.cat', format='ascii')
    #catalog = gaia2file(ra=half_ra, dec=half_dec, size=20., mag_limit=23., output='gaia.cat')
    catalog = gaia3cat(ra=half_ra, dec=half_dec, radius_arcmin=10)
    catalog['raMean'], catalog['decMean'] = catalog['ra'], catalog['dec']
    #fig, ax2 = plt.subplots(1, 1, figsize=(15, 15))
    #x, y = _refdatas.wcs.all_world2pix(catalog['raMean'], catalog['decMean'], 0.)
    #norm = ImageNormalize(_refdatas.data, PercentileInterval(99.))
    #ax2.imshow(_refdatas.data, norm=norm)
    #ax2.plot(x, y, marker='o', mec='r', mfc='none', ls='none')
    #plt.show()
    _, ref_stars = make_psf(_refdatas, catalog, show=show)
    refdatas = _refdatas.copy()
    ##### Calculate the bkg #####
    #data = scidata0.data.byteswap().newbyteorder()
    try:
        _data = refdatas.data.astype(np.float32)
    except:
        _data = refdatas.data.byteswap().newbyteorder()
    bkg = sep.Background(_data)
    ref_global_bkg = bkg.globalback
    #############################
    refine_wcs(refdatas.wcs, ref_stars, catalog)

    if make_psf_temp:
        print(f'make psf for {template_filename}')
        ref_psf, _ = make_psf(refdatas, catalog, show=False)
    else:
        filename = os.path.join(workdir, template_filename.replace('.fits', '.psf.fits'))
        print(f'read the psf for {template_filename}')
        ref_psf = read_with_datasec(filename)


    refdata = assemble_reference([refdatas], scidata.wcs, scidata.shape, ref_global_bkg=ref_global_bkg)
    #_template_filename = os.path.join(workdir, science_filename.replace('.fits', '.temp.fits'))
    #refdata0.write(_template_filename, overwrite=True)
    _science_filename = os.path.join(workdir, 'science_wcs.fits')
    scidata.write(_science_filename, overwrite=True)

    _template_filename = os.path.join(workdir, science_filename.replace('.fits', '.temp.fits'))
    refdata.write(_template_filename, overwrite=True)

    if show:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 15))
        x, y = scidata.wcs.all_world2pix(catalog['raMean'], catalog['decMean'], 0.)

        vmin, vmax = np.percentile(scidata.data, (15., 99.5))
        ax1.imshow(scidata.data, vmin=vmin, vmax=vmax)
        ax1.plot(x, y, marker='o', mec='r', mfc='none', ls='none')

        norm = ImageNormalize(refdata0.data, PercentileInterval(99.))
        ax2.imshow(refdata0.data, norm=norm)
        ax2.plot(x, y, marker='o', mec='r', mfc='none', ls='none')



    # # Subtract the images and view the result

    output_filename = os.path.join(workdir, science_filename.replace('.fits', '.diff.fits'))
    science = ImageClass(scidata, sci_psf.data, saturation=65535)

    refdata.mask[np.isnan(refdata.data)] = True
    #refdata.mask[refdata.data == 0] = True
    reference = ImageClass(refdata, ref_psf.data, refdata.mask, saturation=65535)#, refdata.mask
    #reference.data = np.nan_to_num(reference.data, 
    #        nan=np.nanmedian(np.ma.array(reference.data, mask=reference.mask)))
    #reference.data = np.nan_to_num(reference.data, nan=ref_global_bkg)
    #reference = np.nan_to_num(reference, nan=ref_global_bkg)
    np.nan_to_num(reference.data, nan=ref_global_bkg, copy=False, out=reference.data)
    

    difference = calculate_difference_image(science, reference, show=show, max_iterations=3, sigma_cut = sigma_cut)
    difference_zero_point = calculate_difference_image_zero_point(science, reference)
    normalized_difference = normalize_difference_image(difference, difference_zero_point, science, reference, 'i')
    save_difference_image_to_file(normalized_difference, science, 'i', output_filename)

    #difference_psf = calculate_difference_psf(science, reference, difference_zero_point)
    #save_difference_psf_to_file(difference_psf, output_filename.replace('.fits', '.psf.fits'))

    if show:
        vmin, vmax = np.percentile(difference, (15., 99.5))
        plt.figure(figsize=(7., 15.))
        plt.imshow(difference, vmin=vmin, vmax=vmax)

    # Memory cleanup after each iteration
    import gc
    del scidata0, scidata, sci_psf, refdata, ref_psf, difference, difference_zero_point, normalized_difference
    del catalog# #stars, objects, data_sub, bkg
    gc.collect()
    print(f"Memory cleaned up after processing {science_filename}")