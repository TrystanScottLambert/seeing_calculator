"""
Testing module for extraction seeing profiles.
"""

from typing import Tuple
import unittest
import numpy as np
from astropy.table import Table
from photutils.datasets import make_gaussian_sources_image
import extract_seeing_profile


def create_test_cutout(std: float) -> Tuple[np.ndarray, float]:
    """
    Creates a false star cutout for testing.
    Only free parameter is the std since this is the 
    only important parameter that needs to be recovered.
    """
    sources = Table()
    sources['flux'] = [10000]
    sources['x_mean'] = [30]
    sources['y_mean'] = [30]
    sources['x_stddev'] = [std]
    sources['y_stddev'] = sources['x_stddev']
    sources['theta'] = [0]
    tshape = (60, 60)
    image = (make_gaussian_sources_image(tshape, sources)
)
    return image, std


class TestExtractSeeingProfile(unittest.TestCase):
    """
    Testing class for the extract_seeing_profile module.
    """

    def test_fit2d_gauss(self):
        """
        Testing if the 2d gaussian function is working and recovering several different 
        fake stellar sources.
        """
        test_cutouts = [create_test_cutout(std) for std in np.arange(2, 10, 0.1)]
        for test_cutout in test_cutouts:
            fit = extract_seeing_profile.fit_2d_gauss(test_cutout[0])
            self.assertAlmostEqual(fit.x_stddev.value, test_cutout[1])

if __name__ == '__main__':
    unittest.main()
