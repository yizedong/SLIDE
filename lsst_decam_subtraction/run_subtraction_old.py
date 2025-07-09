#!/usr/bin/env python3
"""
Batch runner for LSST DECam image subtraction.
"""

import os
import subprocess
from astropy.io import fits

def run_batch_subtraction(workdir='./', science_filenames=None, template_filename=None, 
                         ra=58.335054, dec=-48.750303, make_psf_temp=False):
    """
    Run image subtraction on multiple science images.
    
    Parameters
    ----------
    workdir : str
        Working directory
    science_filenames : list
        List of science image filenames
    template_filename : str
        Template image filename (if not downloading)
    ra : float
        Right ascension for reference image query
    dec : float
        Declination for reference image query
    make_psf_temp : bool
        Whether to make PSF for template
    """
    
    if science_filenames is None:
        science_filenames = ['dp1_r_2024120800420.fits']
    
    for science_filename in science_filenames:
        print(f"Processing {science_filename}...")
        
        # Get filter from image header
        try:
            hdul = fits.open(os.path.join(workdir, science_filename))
            header = hdul[0].header
            filt = header.get('FILTBAND', 'r')
            hdul.close()
        except Exception as e:
            print(f"Warning: Could not read filter from header for {science_filename}: {e}")
            filt = 'r'
        
        # Build command
        cmd = [
            'python', '-m', 'lsst_decam_subtraction.main',
            '--targimg', science_filename,
            '--imgdir', workdir,
            '--ra', str(ra),
            '--dec', str(dec),
            '--filter', filt,
            '--download_DES_temp'
        ]
        
        if template_filename:
            cmd.extend(['--tempimg', template_filename])
        
        if make_psf_temp:
            cmd.append('--make_psf_temp')
        
        # Execute command
        print(f"Running: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"Successfully processed {science_filename}")
            if result.stdout:
                print("Output:", result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {science_filename}: {e}")
            if e.stderr:
                print("Error output:", e.stderr)
            continue

if __name__ == "__main__":
    # Configuration
    workdir = './'
    science_filenames = ['dp1_r_2024120800420.fits']
    template_filename = None  # Set to template filename if not downloading
    ra = 58.335054 
    dec = -48.750303
    make_psf_temp = False
    
    # Run batch processing
    run_batch_subtraction(
        workdir=workdir,
        science_filenames=science_filenames,
        template_filename=template_filename,
        ra=ra,
        dec=dec,
        make_psf_temp=make_psf_temp
    )