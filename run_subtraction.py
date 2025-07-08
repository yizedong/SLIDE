import os
from astropy.io import fits

workdir = './'

template_filename = 'dp1_z_2024111700351.fits'
science_filenames = [
                    'dp1_z_2024111800105.fits',
                   'dp1_z_2024111700357.fits',
                    'dp1_z_2024111700356.fits',
                    'dp1_z_2024120800434.fits',
                    'dp1_z_2024120800435.fits'
                   ]
template_filename = 'dp1_r_2024111700344.fits'
science_filenames = [#'dp1_r_2024112600180.fits',
                    'dp1_r_2024120200090.fits',
                    'dp1_r_2024111900342.fits',
                    'dp1_r_2024111800097.fits',
                    'dp1_r_2024111700369.fits',
                    'dp1_r_2024112600182.fits',
                    'dp1_r_2024120800420.fits',
                    'dp1_r_2024112600183.fits',
                    'dp1_r_2024111700368.fits'
                    ]

science_filenames = ['dp1_r_2024112600180.fits',
                        'dp1_r_2024120200090.fits',
                        'dp1_r_2024111900342.fits',
                        'dp1_r_2024111800097.fits',
                        'dp1_r_2024111700369.fits',
                        'dp1_r_2024112600182.fits',
                        'dp1_r_2024120800420.fits',
                        'dp1_r_2024111700344.fits',
                        'dp1_r_2024112600183.fits',
                        'dp1_r_2024111700368.fits',
                        'dp1_g_2024120800422.fits', 
                        'dp1_g_2024120800423.fits',
                        'dp1_i_2024112500297.fits', 
                        'dp1_i_2024112600167.fits', 
                        'dp1_i_2024111700360.fits',
                        'dp1_z_2024111800105.fits',
                        'dp1_z_2024111700357.fits',
                        'dp1_z_2024111700356.fits',
                        'dp1_z_2024111700351.fits',
                        'dp1_z_2024120800434.fits',
                        'dp1_z_2024120800435.fits'
                        ]
#template_filename = 'dp1_i_2024111700360.fits'
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