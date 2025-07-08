import os
from astropy.io import fits

workdir = './'

science_filenames = ['dp1_r_2024120800420.fits']
make_psf_temp = False
SN_ra = 58.335054 
SN_dec = -48.750303

for science_filename in science_filenames:
    hdul = fits.open(workdir + science_filename)
    header = hdul[0].header
    filter = header['FILTBAND']
    command = f'python decam_image_subtraction.py --targimg {science_filename}\
        --tempimg {template_filename} --imgdir {workdir} --ra {SN_ra} --dec {SN_dec} --filter {filter} --download_DES_temp --make_psf_temp' #--make_psf_temp
    print(command)
    os.system(command)