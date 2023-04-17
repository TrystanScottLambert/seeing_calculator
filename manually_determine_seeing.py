"""
Program which allows the user to 
manually determine the seeing of the image by manually
selecting bright, but non-saturated stars.
"""

import numpy as np
import pylab as plt
from astropy.io import fits
from astropy.visualization import ZScaleInterval

from extract_seeing_profile import fit_2d_gauss, get_average_fwhm

BOX_WIDTH = 15 # pixels

def stamp_cut_out(data: np.ndarray, x_pos:float, y_pos:float) -> np.ndarray:
    """
    Makes a standard stamp cut out of an image, using a 15x15 pixel cut
    around eaach center
    """
    cut_out = data[
        int(y_pos)-BOX_WIDTH:int(y_pos)+BOX_WIDTH, int(x_pos)-BOX_WIDTH:int(x_pos)+BOX_WIDTH]
    return cut_out



class FitsPlot:
    """Class for the main image plotting."""
    def __init__(self, infile: str) -> None:
        """Initializing."""
        hdul = fits.open(infile)
        self.data = hdul[0].data
        z_interval = ZScaleInterval()
        self.vmin, self.vmax = z_interval.get_limits(self.data)
        self.cut_outs = []
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
    
    def plot_whole(self):
        """Plots the entire fits file."""

        self.ax.imshow(self.data, vmin=self.vmin, vmax=self.vmax)
        _ = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        plt.show()
    
    def show_cut_out(self, data: np.ndarray) -> None:
        """Plots what the cutout looks like"""
        side_fig = plt.figure()
        side_ax = side_fig.add_subplot(111)
        side_ax.imshow(data)
        cid = side_fig.canvas.mpl_connect('button_press_event', self.onclick_accept)
        plt.show()

    def onclick(self, event):
        """
        Interacting with the main plot. On dbl click ...
        """
        if event.dblclick:
            ix, iy = event.xdata, event.ydata
            first_cut_out = stamp_cut_out(self.data, ix, iy)
            model = fit_2d_gauss(first_cut_out)
            x_model, y_model = model.x_mean.value, model.y_mean.value

            offset_x = x_model-BOX_WIDTH
            offset_y = y_model-BOX_WIDTH
            center_x, center_y = ix+offset_x,  iy+offset_y
            self.ax.scatter(center_x, center_y, color='r')
            self.cut_out = stamp_cut_out(self.data, center_x, center_y)
            self.show_cut_out(self.cut_out)

    def onclick_accept(self, event) -> None:
        """
        Clicking to accept the star
        """
        if event.dblclick:
            print('DBL CLICK')
            self.cut_outs.append(self.cut_out)

    def calculate_seeing(self):
        """
        Determines the average fwhm for every cutout image that was saved.
        """
        if len(self.cut_outs) == 0:
            raise ValueError("No cutouts have been selected.")

        avg_fwhms = [get_average_fwhm(cut_out) for cut_out in self.cut_outs]
        return np.mean(avg_fwhms)


if __name__ == '__main__':
    INFILE = '../DECAM_analysis/correct_stacks/N964/n964.fits'
    test = FitsPlot(INFILE)
    test.plot_whole()
    seeing = test.calculate_seeing()
