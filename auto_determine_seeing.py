"""
Module to automatically determine the seeing of a fits image.
"""

from typing import List
import numpy as np
from pylab import plt
from astropy.io import fits
from source_identification import find_appropriate_sources
from extract_seeing_profile import get_average_fwhm

def stamp_cut_out(data: np.ndarray, x_pos:float, y_pos:float) -> np.ndarray:
    """
    Makes a standard stamp cut out of an image, using a 30x30 pixel cut
    around eaach center
    """
    cut_out = data[int(y_pos)-30:int(y_pos)+30, int(x_pos)-30:int(x_pos)+30]

    return cut_out

def get_avg_fwhm_of_positions(data: np.ndarray, x_positions: List[float], y_positions: List[float]):
    """
    Takes a list of centers and determines fwhm of the stars at each of the positions.
    """
    fwhms = []
    for i, x_pos in enumerate(x_positions):
        cut_out = stamp_cut_out(data, x_pos, y_positions[i])
        if cut_out.shape[0] == cut_out.shape[1]: # Ignore case were one axis is 0 (edge cases)
            print("WE DID IT")
            plt.imshow(cut_out)
            plt.show()
            fwhm = get_average_fwhm(cut_out)
            fwhms.append(fwhm)
    return np.mean(fwhm)

def determine_seeing(data: np.ndarray):
    """
    Main function for automatically determining the seeing of a fits image.
    Note that the seeing will be returned as pixels. Multiply by the pixel scale
    to get the real-world values.
    """
    sources = find_appropriate_sources(data)
    x_positions = np.array(sources['xcentroid'].value).astype(int)
    y_positions = np.array(sources['ycentroid'].value).astype(int)
    fwhm = get_avg_fwhm_of_positions(data, x_positions, y_positions)
    return fwhm

if __name__ == '__main__':
    INFILE = '../DECAM_analysis/correct_stacks/N964/n964.fits'
    hdu = fits.open(INFILE)
    seeing = determine_seeing(hdu[0].data)
    print(seeing)
