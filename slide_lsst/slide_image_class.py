# This is a modified version of the ImageClass from https://github.com/dguevel/PyZOGY/blob/master/PyZOGY/image_class.py.

import numpy as np
from astropy.io import fits
from PyZOGY import util
from scipy.ndimage import gaussian_filter

def fast_interpolate_bad_pixels(image, mask, sigma=3):
    """Faster version using Gaussian blur only on masked pixels."""
    blurred = gaussian_filter(image, sigma=sigma)
    filled_image = image.copy()
    filled_image[mask] = blurred[mask]
    return filled_image

class FASTImageClass(np.ndarray):
    """
    Handles images used by subtraction functions.

    This class handles image data, PSF data, mask data, and other parameters needed for 
    subtraction. It is a subclass of numpy ndarray, so the image data can be accessed
    like a normal numpy array.
    
    Args:
        image_filename (str): Name of the FITS file containing the image.
        psf_filename (str, optional): Name of the FITS file containing the PSF.
        mask_filename (str, optional): Name of the FITS file containing the bad pixel 
        mask array with 1 indicating masking and 0 indicating no masking.
        n_stamps (int, optional): Number of stamps to use for background estimation.
        saturation (float, optional): Maximum usable value in the FITS image.
        read_noise (float, optional): Read noise of the FITs image.

    Attributes:
        header (astropy.io.fits.Header): Header from the image FITS file.
        raw_image (numpy.ndarray): Unaltered data from the image FITS file.
        raw_psf (numpy.ndarray): Unaltered data frin tge PSF FITS file.
        background_std (float): Standard deviation of the image.
        image_filename (str): Filename of the image FITS.
        psf_filename (str): Filename of the PSF FITS.
        saturation (float): Maximum usable value in the FITS image.
        mask (numpy.ndarray): Bad pixel mask for the image.
        psf (numpy.ndarray): PSF after shifting, resizing, and normalization.
        zero_point (float): Flux based zero point of the image
    """

    def __new__(cls, image_filename, psf_filename, mask_filename=None, n_stamps=1,
                saturation=np.inf, variance=None, read_noise=0, registration_noise=(0, 0)):
        if isinstance(image_filename, str):  # filename
            raw_image, header = fits.getdata(image_filename, header=True)
        elif hasattr(image_filename, 'header'):  # FITS HDU or astropy CCDData
            raw_image = image_filename.data
            header = image_filename.header
            image_filename = 'IMAGE_IN_MEMORY'
        else:  # numpy array
            raw_image = image_filename
            header = None
            image_filename = 'IMAGE_IN_MEMORY'
        if isinstance(psf_filename, str):
            raw_psf = fits.getdata(psf_filename)
        else:
            raw_psf = psf_filename
            psf_filename = 'PSF_IN_MEMORY'
        psf = util.center_psf(util.resize_psf(raw_psf, raw_image.shape), fname=image_filename) / np.sum(raw_psf)
        if isinstance(mask_filename, str):
            mask = fits.getdata(mask_filename)
        else:
            mask = mask_filename
        mask = util.mask_saturated_pix(raw_image, saturation, mask, image_filename)
        masked_image = np.ma.array(raw_image, mask=mask)
        background_std, background_counts = util.fit_noise(masked_image, n_stamps=n_stamps, fname=image_filename)
        #image_data = util.interpolate_bad_pixels(masked_image, fname=image_filename) - background_counts
        image_data = fast_interpolate_bad_pixels(masked_image.data, masked_image.mask, sigma=3)

        obj = np.asarray(image_data).view(cls)
        obj.header = header
        obj.raw_image = raw_image
        obj.raw_psf = raw_psf
        obj.background_std = background_std
        obj.background_counts = background_counts
        obj.image_filename = image_filename
        obj.psf_filename = psf_filename
        obj.saturation = saturation
        obj.mask = mask
        obj.psf = psf
        obj.zero_point = 1.
        obj.variance = variance
        obj.read_noise = read_noise
        obj.registration_noise = registration_noise

        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.raw_image = getattr(obj, 'raw_image', None)
        self.header = getattr(obj, 'header', None)
        self.raw_psf = getattr(obj, 'raw_psf', None)
        self.background_std = getattr(obj, 'background_std', None)
        self.background_counts = getattr(obj, 'background_counts', None)
        self.image_filename = getattr(obj, 'image_filename', None)
        self.psf_filename = getattr(obj, 'psf_filename', None)
        self.saturation = getattr(obj, 'saturation', None)
        self.mask = getattr(obj, 'mask', None)
        self.psf = getattr(obj, 'psf', None)
        self.zero_point = getattr(obj, 'zero_point', None)
        self.variance = getattr(obj, 'variance', None)
        self.read_noise = getattr(obj, 'read_noise', None)
        self.registration_noise = getattr(obj, 'registration_noise', None)