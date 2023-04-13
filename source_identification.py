"""
Module for identiyfing non saturated sources in the image. 
"""

import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder

def find_sources(data: np.ndarray) -> Table:
    """Identifies the bright sources in an image."""
    _, median, std = sigma_clipped_stats(data, sigma=3.0)
    daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
    sources = daofind(data - median)
    return sources

def select_non_saturated_sources(source_table: Table):
    """
    Determines the sources that are not saturated and are
    not too dim in order to determine the fwhm. To do this
    we take the values between median - median/4 to median 
    + median/4.
    """
    fluxes = source_table['flux'].value
    median = np.median(fluxes)
    cut = np.where((fluxes > median - median/4) & (fluxes < median + median/4))[0]
    return cut

def find_appropriate_sources(data: np.ndarray) -> Table:
    """
    Returns only good sources in a given image.
    """
    all_sources = find_sources(data)
    idx_good = select_non_saturated_sources(all_sources)
    return all_sources[idx_good]

if __name__ == '__main__':
    INFILE = '../DECAM_analysis/correct_stacks/N964/n964.fits'
    hdu = fits.open(INFILE)
    good_sources = find_appropriate_sources(hdu[0].data)
