"""
Module to perform aperture photometry on a single stellar source. 
I.e., place several annuli over a 2d source and extract a profile. 
"""

import numpy as np
from astropy.modeling import models, fitting

def fit_2d_gauss(data: np.ndarray) -> models.Gaussian2D:
    """
    Models a 2 Gaussian and fits it to the cut out source.
    The source has to be a stamp cut out. Function cannot be
    used on a full astronomy image.
    """
    fit = fitting.LevMarLSQFitter()
    # Initial Guesses based on the data
    y_0, x_0 = np.unravel_index(np.argmax(data), data.shape)
    sigma = np.std(data)
    amp = np.max(data)

    model = models.Gaussian2D(amp, x_0, y_0, sigma, sigma)

    y_i, x_i = np.indices(data.shape)

    fitted_model = fit(model, x_i, y_i, data)
    return fitted_model

def get_average_fwhm(data: np.ndarray):
    """
    Determines the average standard deviation of the fitted 2D Gaussian
    and then converts this into FWHM
    """
    gaussian_2d = fit_2d_gauss(data)
    avg_fwhm = np.mean([gaussian_2d.x_fwhm, gaussian_2d.y_fwhm])
    return avg_fwhm
