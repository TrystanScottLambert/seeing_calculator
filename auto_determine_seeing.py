"""
Automatically estimating the seeing.
"""

import numpy as np
from astropy.io import fits
from astropy.stats import gaussian_fwhm_to_sigma

from source_identification import find_appropriate_sources
from extract_seeing_profile import get_average_fwhm


BOX_WIDTH = 15 # pixels

def stamp_cut_out(data: np.ndarray, x_pos:float, y_pos:float) -> np.ndarray:
    """
    Makes a standard stamp cut out of an image, using a 15x15 pixel cut
    around eaach center
    """
    cut_out = data[
        int(y_pos)-BOX_WIDTH:int(y_pos)+BOX_WIDTH, int(x_pos)-BOX_WIDTH:int(x_pos)+BOX_WIDTH]
    return cut_out

def calc_seeing(data: np.ndarray, fwhm:float, std_above_background:float) -> list[float]:
    """Returns a list of fwhm measurements for in pixels."""
    sources = find_appropriate_sources(data, fwhm, std_above_background)
    x_positions = sources['xcentroid'].value
    y_positions = sources['ycentroid'].value
    fwhms = []
    for i, _ in enumerate(x_positions):
        try:
            postage_stamp = stamp_cut_out(data, x_positions[i], y_positions[i])
            fwhms.append(get_average_fwhm(postage_stamp))
        except:
            pass
    return np.mean(fwhms) * gaussian_fwhm_to_sigma


if __name__ == '__main__':
    INFILE = '../DECAM_analysis/correct_stacks/N964/n964.fits'
    #INFILE = '/home/tlambert/Desktop/IMACS_photometry/imacs_data/night_1_theli.fits'
    hdu = fits.open(INFILE)
    decam_n_seeing = calc_seeing(hdu[0].data, fwhm=7.5, std_above_background=25)
    print(decam_n_seeing)
