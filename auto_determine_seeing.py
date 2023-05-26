"""
Automatically estimating the seeing.
"""

import numpy as np
from astropy.io import fits
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils.background import Background2D, MedianBackground
from source_identification import find_appropriate_sources
from extract_seeing_profile import get_average_fwhm

sigma_clip = SigmaClip(sigma=3.0)
bkg_estimator = MedianBackground()


BOX_WIDTH = 15 # pixels

def stamp_cut_out(data: np.ndarray, x_pos:float, y_pos:float) -> np.ndarray:
    """
    Makes a standard stamp cut out of an image, using a 15x15 pixel cut
    around eaach center
    """
    cut_out = data[
        int(y_pos)-BOX_WIDTH:int(y_pos)+BOX_WIDTH, int(x_pos)-BOX_WIDTH:int(x_pos)+BOX_WIDTH]
    return cut_out

def calc_seeing(data: np.ndarray, fwhm:float, min_val:float, plot: bool = False) -> float:
    """Returns a list of fwhm measurements for in pixels."""
    sources = find_appropriate_sources(data, fwhm, min_val, plot)
    x_positions = sources['xcentroid'].value
    y_positions = sources['ycentroid'].value
    fwhms = []
    for i, _ in enumerate(x_positions):
        try:
            postage_stamp = stamp_cut_out(data, x_positions[i], y_positions[i])
            fwhms.append(get_average_fwhm(postage_stamp))
        except:
            pass

    fwhms = np.array(fwhms)
    fwhms = fwhms[(fwhms>0.5) & (fwhms<20)]  # More than 20 pixels is just unphysical
    return sigma_clipped_stats(fwhms)[1]


if __name__ == '__main__':
    import pylab as plt
    from astropy.visualization import ZScaleInterval
    INFILES = [
    #'../IMACS_photometry/imacs_data/night_1_theli.fits',
    #'../IMACS_photometry/imacs_data/night_2_theli.fits',
    '../DECAM_analysis/correct_stacks/N964/z.fits',
    '../DECAM_analysis/correct_stacks/N964/.fits',
    '../DECAM_analysis/correct_stacks/N964/z.fits',
    ]

    for file in INFILES:
        hdu = fits.open(file)
        len_y, len_x = hdu[0].data.shape
        center_y = int(len_y/2)
        center_x = int(len_x/2)
        zscale = ZScaleInterval()
        v_min, v_max = zscale.get_limits(hdu[0].data)
        PAD = 6000
        data = hdu[0].data[center_y-PAD: center_y+PAD, center_x-PAD: center_x+PAD]
        bkg = Background2D(data, data.shape, filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        data -= bkg.background
        plt.imshow(data)
        zscale = ZScaleInterval()
        v_min, v_max = zscale.get_limits(data)
        seeing = calc_seeing(data, fwhm=5, min_val=10)
        print(file.split('/')[-1], seeing)
