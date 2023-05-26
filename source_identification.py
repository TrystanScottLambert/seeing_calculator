"""
Module for identiyfing non saturated sources in the image. 
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval
from photutils.detection import IRAFStarFinder
from photutils.aperture import CircularAperture


def find_center(image: np.ndarray, factor: float) -> np.ndarray:
    """
    Returns a mask representing the 1./factor central region of the image.
    """

    height, width = image.shape
    center_x = int(image.shape[1]/2)
    center_y = int(image.shape[0]/2)
    box_height = int(height/factor)
    box_width = int(width/factor)

    mask = np.zeros(image.shape, dtype=bool)
    mask[:,:] = True
    mask[center_y-box_height:center_y + box_height,center_x-box_width: center_x+box_width] = 0.
    return mask

def find_sources(
        data: np.ndarray, fwhm: float, min_val: float, central_factor: float = 4.
        ) -> Table:
    """Identifies the bright sources in an image."""
    _, median, _ = sigma_clipped_stats(data)
    daofind = IRAFStarFinder(fwhm=fwhm, threshold=min_val)
    if central_factor:
        msk = find_center(data, central_factor)
        sources = daofind(data - median, mask=msk)
    else:
        sources = daofind(data - median)
    return sources

def select_non_saturated_sources(source_table: Table):
    """
    Determines the sources that are not saturated and are
    not too dim in order to determine the fwhm. 
    """
    fluxes = source_table['flux'].value
    median = np.median(fluxes)
    cut = np.where(fluxes < 2*median)[0]
    return cut

def find_appropriate_sources(
        data: np.ndarray, fwhm: float, min_val: float,
        central_factor: float = 4, plot=False
        ) -> Table:
    """
    Returns only good sources in a given image.
    """
    all_sources = find_sources(data, fwhm, min_val, central_factor=central_factor)
    idx_good = select_non_saturated_sources(all_sources)
    good_sources = all_sources[idx_good]

    if plot:
        positions = np.transpose((good_sources['xcentroid'], good_sources['ycentroid']))
        apertures = CircularAperture(positions, r=4.0)
        zscale = ZScaleInterval()
        v_min, v_max = zscale.get_limits(data)
        plt.imshow(data, cmap='Greys', vmin=v_min, vmax=v_max)
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        plt.show()
    return good_sources
