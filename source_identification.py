"""
Module for identiyfing non saturated sources in the image. 
"""

import numpy as np
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from photutils.detection import IRAFStarFinder

def find_sources(data: np.ndarray, fwhm: float, std_above_background: float) -> Table:
    """Identifies the bright sources in an image."""
    _, median, std = sigma_clipped_stats(data)
    daofind = IRAFStarFinder(fwhm=fwhm, threshold=std_above_background*std)
    sources = daofind(data - median)
    return sources

def select_non_saturated_sources(source_table: Table):
    """
    Determines the sources that are not saturated and are
    not too dim in order to determine the fwhm. To do this
    we take the values between median + median/2 to median.
    """
    fluxes = source_table['flux'].value
    median = np.median(fluxes)
    cut = np.where(fluxes < 2*median)[0]
    return cut

def find_appropriate_sources(data: np.ndarray, fwhm: float, std_above_background: float) -> Table:
    """
    Returns only good sources in a given image.
    """
    all_sources = find_sources(data, fwhm, std_above_background)
    idx_good = select_non_saturated_sources(all_sources)
    return all_sources[idx_good]
